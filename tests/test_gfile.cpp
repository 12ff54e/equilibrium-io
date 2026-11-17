#include "equilibrium/gfile.hpp"
#include "equilibrium/record_validation.hpp"
#include "equilibrium/writer_helpers.hpp"
#include "gfile_fixture.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double actual,
           double expected,
           double tolerance,
           const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr << message << ": " << actual << " expected " << expected
                  << '\n';
        throw std::runtime_error(message);
    }
}
void rejects(const auto& function, const char* message) {
    try {
        function();
    } catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}
equilibrium::GfileInput convert(
    const gfile_test::Settings& source = {},
    const equilibrium::GfileInputSettings& settings = {}) {
    std::istringstream stream(gfile_test::fixture(source));
    return equilibrium::gfile_to_cumes_input(stream, settings);
}
void shapes() {
    using equilibrium::LcfsPoint;
    const std::vector<LcfsPoint> points{
        {4, .2}, {3.1, 2}, {2.7, 2}, {2, .2}, {3.3, -1.5}};
    auto verify = [&](const auto& polygon, double shift_r, double shift_z) {
        const auto s = equilibrium::calculate_lcfs_shape(polygon, .2 + shift_z);
        close(s.r_center, 3 + shift_r, 1e-14, "shape radial center");
        close(s.minor_radius, 1, 1e-14, "shape minor radius");
        close(s.elongation_upper, 1.8, 1e-14, "upper elongation");
        close(s.elongation_lower, 1.7, 1e-14, "lower elongation");
        close(s.triangularity_upper, .1, 1e-14, "flat upper triangularity");
        close(s.triangularity_lower, -.3, 1e-14,
              "negative lower triangularity");
    };
    verify(points, 0, 0);
    auto duplicate = points;
    duplicate.insert(duplicate.begin() + 1, duplicate.front());
    duplicate.push_back(duplicate.front());
    verify(duplicate, 0, 0);
    duplicate.insert(duplicate.begin() + 3, {3.0, 2});
    verify(duplicate, 0, 0);
    std::reverse(duplicate.begin(), duplicate.end());
    verify(duplicate, 0, 0);
    for (auto& p : duplicate) {
        p.r += 10;
        p.z += 7;
    }
    verify(duplicate, 10, 7);
    auto reflected = points;
    for (auto& p : reflected) p.z = .4 - p.z;
    const auto reflection = equilibrium::calculate_lcfs_shape(reflected, .2);
    close(reflection.elongation_upper, 1.7, 1e-14, "reflected elongation");
    close(reflection.triangularity_upper, -.3, 1e-14,
          "reflected triangularity");
    rejects([&] { equilibrium::calculate_lcfs_shape(points, 2); },
            "invalid shape reference accepted");
    rejects(
        [&] {
            equilibrium::calculate_lcfs_shape(
                std::vector<LcfsPoint>{{3, 0}, {3, 1}, {3, -1}}, 0);
        },
        "zero width accepted");
}
void conventions(const std::filesystem::path& directory) {
    using namespace gfile_test;
    for (int offset : {0, 10})
        for (int index = 1; index <= 8; ++index)
            for (int ip : {-1, 1})
                for (int bt : {-1, 1}) {
                    Settings source;
                    source.cocos = index + offset;
                    source.current_sign = ip;
                    source.field_sign = bt;
                    const auto converted = convert(source);
                    const auto& p = converted.input;
                    check(converted.lcfs.cocos == source.cocos,
                          "COCOS resolution mismatch");
                    check(bool(equilibrium::validate_input_record(p)),
                          "invalid normalized input");
                    check(p.lasym && p.ntor == 0 && p.ncurr == 0 &&
                              p.nzeta == 1 && p.ntheta == 30,
                          "incorrect input defaults");
                    close(p.phiedge,
                          bt * 2 * std::numbers::pi * PSI * (Q0 + .5 * Q1),
                          2e-8, "signed full toroidal flux");
                    close(p.curtor, ip * 2.e6, 1e-6, "signed physical current");
                    close(p.zaxis_c[0], source.axis_z, 1e-12,
                          "absolute magnetic axis Z");
                    for (std::size_t k = 0; k < p.am.size(); ++k) {
                        const double s =
                            static_cast<double>(k) / (p.am.size() - 1);
                        const double u = poloidal_fraction(s);
                        close(p.am[k], P0 * (1 - .5 * u), 1.e-4,
                              "pressure remap to toroidal flux");
                        close(p.ai[k], -ip * bt / (Q0 + Q1 * u), 2.e-9,
                              "signed iota remap");
                    }
                    const auto shape =
                        equilibrium::calculate_lcfs_shape(converted.lcfs);
                    close(shape.elongation_upper, B / A, 1e-8,
                          "gfile upper elongation");
                    close(shape.elongation_lower, B / A, 1e-8,
                          "gfile lower elongation");
                    close(shape.triangularity_upper, 0, 1e-8,
                          "ellipse triangularity");
                    check(converted.boundary_max_error < 1.e-4 &&
                              converted.poloidal_flux_relative_error < 1.e-4,
                          "conversion accuracy");
                    if (source.cocos == 5) {
                        const auto name = std::string("analytic-") +
                                          (ip > 0 ? "ip+" : "ip-") +
                                          (bt > 0 ? "bt+" : "bt-");
                        std::ofstream(directory / (name + ".g"))
                            << fixture(source);
                        const auto status = equilibrium::write_gfile_input(
                            (directory / (name + ".json")).string(), converted);
                        check(bool(status), "writing generated JSON failed");
                    }
                }
    Settings symmetric;
    symmetric.axis_z = 0;
    const auto converted = convert(symmetric);
    check(!converted.input.lasym && converted.input.zaxis_c.empty(),
          "symmetric input marked asymmetric");
    std::ofstream(directory / "symmetric.g") << fixture(symmetric);
    check(bool(equilibrium::write_gfile_input(
              (directory / "symmetric.json").string(), converted)),
          "write symmetric JSON");
    Settings alternate;
    alternate.close_boundary = true;
    alternate.reverse_boundary = true;
    alternate.fortran_d = true;
    const auto other = convert(alternate);
    close(other.input.phiedge, convert().input.phiedge, 1e-14,
          "Fortran D/clockwise boundary changes flux");
    Settings unlabeled;
    unlabeled.cocos = 11;
    unlabeled.header = "EFIT unlabelled full flux";
    equilibrium::GfileInputSettings forced;
    forced.cocos = 11;
    close(convert(unlabeled, forced).input.phiedge, convert().input.phiedge,
          1e-8, "explicit COCOS override");
}
double interpolate(const std::vector<double>& values, double s) {
    const double x = s * (values.size() - 1);
    const auto k = std::min(static_cast<std::size_t>(x), values.size() - 2);
    return std::lerp(values[k], values[k + 1], x - static_cast<double>(k));
}
void refinement_and_cropping() {
    using namespace gfile_test;
    equilibrium::GfileInputSettings coarse;
    coarse.profile_points = 9;
    coarse.poloidal_flux_tolerance = 1;
    auto fine = coarse;
    fine.profile_points = 65;
    const auto a = convert({}, coarse), b = convert({}, fine);
    const auto error = [](const auto& result) {
        double maximum = 0;
        for (int k = 0; k <= 1000; ++k) {
            const double s = k / 1000.0, u = poloidal_fraction(s);
            maximum = std::max(
                maximum,
                std::abs(interpolate(result.input.ai, s) + 1 / (Q0 + Q1 * u)));
            maximum = std::max(
                maximum,
                std::abs(interpolate(result.input.am, s) / P0 - (1 - .5 * u)));
        }
        return maximum;
    };
    check(error(b) < error(a) / 20, "profile refinement did not converge");
    check(b.poloidal_flux_relative_error < a.poloidal_flux_relative_error / 20,
          "poloidal-flux error did not converge");
    coarse.mpol = 4;
    coarse.boundary_tolerance = 1;
    const auto rough = convert({}, coarse);
    check(b.boundary_max_error < rough.boundary_max_error / 10,
          "boundary Fourier refinement did not converge");
    equilibrium::GfileInputSettings crop;
    crop.psi_ratio = .64;
    const auto cropped = convert({}, crop);
    close(cropped.input.phiedge,
          2 * std::numbers::pi * PSI * (Q0 * .64 + .5 * Q1 * .64 * .64), 1e-8,
          "cropped toroidal flux");
    close(cropped.input.am.back(), P0 * (1 - .5 * .64), 1e-5,
          "cropped edge pressure");
    double rout = 0;
    for (double coefficient : cropped.input.rbc_value) rout += coefficient;
    close(rout, R0 + .8 * A, 1e-5, "cropped ellipse outer radius");
    close(equilibrium::calculate_lcfs_shape(cropped.lcfs).minor_radius, A, 1e-8,
          "shape must use original LCFS");
    for (int index : {2, 11, 16})
        for (int ip : {-1, 1})
            for (int bt : {-1, 1}) {
                Settings source;
                source.cocos = index;
                source.current_sign = ip;
                source.field_sign = bt;
                const auto signed_crop = convert(source, crop);
                close(signed_crop.input.phiedge, bt * cropped.input.phiedge,
                      1e-8, "cropped COCOS flux invariance");
                close(signed_crop.input.ai.back(), -ip * bt / (Q0 + Q1 * .64),
                      1e-9, "cropped COCOS iota sign");
                for (std::size_t m = 0; m < cropped.input.rbcc.size(); ++m)
                    close(signed_crop.input.rbcc[m], cropped.input.rbcc[m],
                          1e-8, "cropped COCOS geometry invariance");
            }
    Settings tilted;
    tilted.tilt = .3;
    const auto asymmetric = convert(tilted);
    check(
        std::any_of(asymmetric.input.rbsc.begin(), asymmetric.input.rbsc.end(),
                    [](double v) { return std::abs(v) > 1e-3; }),
        "tilted boundary lost R sine coefficients");
    for (int k = 0; k < 1024; ++k) {
        const double t = 2 * std::numbers::pi * k / 1024;
        double r = 0, z = 0;
        for (int m = 0; m < asymmetric.input.mpol; ++m) {
            r += asymmetric.input.rbcc[m] * std::cos(m * t) +
                 asymmetric.input.rbsc[m] * std::sin(m * t);
            z += asymmetric.input.zbcc[m] * std::cos(m * t) +
                 asymmetric.input.zbsc[m] * std::sin(m * t);
        }
        r -= R0;
        z -= tilted.axis_z;
        close(r * r / (A * A) + z * z / (B * B) +
                  2 * tilted.tilt * r * z / (A * B),
              1, 2e-4, "tilted analytic LCFS reconstruction");
    }
    Settings unbracketed;
    unbracketed.grid_flux_scale = .5;
    rejects([&] { convert(unbracketed, crop); }, "unbracketed crop accepted");
    auto too_coarse = coarse;
    too_coarse.profile_points = 2;
    too_coarse.poloidal_flux_tolerance = 1e-3;
    rejects([&] { convert({}, too_coarse); },
            "coarse iota flux error accepted");
    coarse.boundary_tolerance = 1e-10;
    rejects([&] { convert({}, coarse); }, "poor boundary fit accepted");
}
void malformed_and_geometry_only(const std::filesystem::path& directory) {
    using namespace gfile_test;
    Settings zero;
    zero.zero_profiles = true;
    std::istringstream geometry(fixture(zero));
    const auto lcfs = equilibrium::read_gfile_lcfs(geometry);
    check(equilibrium::calculate_lcfs_shape(lcfs).elongation_upper > 1,
          "zero-q shape read failed");
    rejects([&] { convert(zero); }, "zero q/F accepted for conversion");
    Settings nonstar;
    nonstar.nonstar = true;
    std::istringstream bent(fixture(nonstar));
    check(equilibrium::read_gfile_lcfs(bent).points.size() == 128,
          "nonstar shape-only read failed");
    rejects([&] { convert(nonstar); },
            "nonstar boundary accepted for conversion");
    Settings reversed_q;
    reversed_q.flip_q = true;
    rejects([&] { convert(reversed_q); }, "q reversal accepted");
    for (const auto& text :
         std::vector<std::string>{"", fixture().substr(0, 180),
                                  fixture().substr(0, fixture().size() - 40)}) {
        rejects(
            [&] {
                std::istringstream in(text);
                equilibrium::read_gfile_lcfs(in);
            },
            "truncated file accepted");
    }
    for (const auto& field : std::vector<std::string>{"  -1", "9999", "abcd"}) {
        auto text = fixture();
        text.replace(52, 4, field);
        rejects(
            [&] {
                std::istringstream in(text);
                equilibrium::read_gfile_lcfs(in);
            },
            "bad grid dimension accepted");
    }
    auto bad_count = fixture();
    const auto pos = bad_count.find("128 0\n");
    check(pos != std::string::npos, "fixture count missing");
    bad_count.replace(pos, 6, "0 0\n");
    rejects(
        [&] {
            std::istringstream in(bad_count);
            equilibrium::read_gfile_lcfs(in);
        },
        "zero boundary count accepted");
    Settings bad_tag;
    bad_tag.header = "COCOS=99";
    rejects([&] { convert(bad_tag); }, "invalid COCOS tag accepted");
    equilibrium::GfileInputSettings invalid;
    invalid.mpol = 17;
    rejects([&] { convert({}, invalid); }, "unsupported cuMES mpol accepted");
    equilibrium::GfileInputSettings limits;
    limits.ntheta = 256;
    limits.stage.ns = 512;
    limits.stage.ftol = 1e-16;
    const auto at_limits = convert({}, limits);
    check(at_limits.input.ntheta == 256 &&
              at_limits.input.stages.front().ns == 512,
          "valid cuMES limits rejected");
    limits.ntheta = 258;
    rejects([&] { convert({}, limits); }, "cuMES ntheta upper limit ignored");
    limits.ntheta = 256;
    limits.stage.ns = 513;
    rejects([&] { convert({}, limits); }, "cuMES ns upper limit ignored");
    limits.stage.ns = 512;
    limits.stage.ftol = 1e-17;
    rejects([&] { convert({}, limits); },
            "cuMES double tolerance floor ignored");
    Settings degenerate;
    degenerate.axis_z = 0;
    degenerate.degenerate_angles = true;
    rejects(
        [&] {
            std::istringstream in(fixture(degenerate));
            equilibrium::read_gfile_lcfs(in);
        },
        "libmeq normalized-angle collapse was not caught");
    const auto valid = convert();
    const auto collision_path = (directory / "collision.json").string();
    std::ofstream(collision_path) << "destination";
    const auto previous = equilibrium::io_detail::temp_path_for(collision_path);
    const auto dot = previous.rfind('.');
    const auto collision =
        previous.substr(0, dot + 1) +
        std::to_string(std::stoul(previous.substr(dot + 1)) + 1);
    std::ofstream(collision) << "other writer";
    check(!equilibrium::write_gfile_input(collision_path, valid),
          "exclusive-create collision accepted");
    std::ifstream collided(collision);
    std::string collision_bytes((std::istreambuf_iterator<char>(collided)), {});
    check(collision_bytes == "other writer",
          "exclusive-create failure deleted an unowned file");
    std::ifstream destination(collision_path);
    std::string destination_bytes;
    destination >> destination_bytes;
    check(destination_bytes == "destination",
          "exclusive-create failure modified destination");
    std::filesystem::remove(collision);
    auto result = convert();
    const auto path = (directory / "preserved.json").string();
    std::ofstream(path) << "original";
    result.input.ai.clear();
    check(!equilibrium::write_gfile_input(path, result),
          "invalid JSON output succeeded");
    std::ifstream saved(path);
    std::string content;
    saved >> content;
    check(content == "original",
          "invalid conversion clobbered existing output");
    rejects([&] { equilibrium::read_gfile_lcfs(directory / "missing.g"); },
            "missing file accepted");
}
}  // namespace
int main(int argc, char** argv) {
    try {
        const auto directory = argc > 1 ? std::filesystem::path(argv[1])
                                        : std::filesystem::current_path();
        shapes();
        conventions(directory);
        refinement_and_cropping();
        malformed_and_geometry_only(directory);
        std::cout << "G-file conversion and LCFS shapes: PASS (64 "
                     "COCOS/polarity cases)\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
