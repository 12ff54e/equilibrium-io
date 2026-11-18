#include "equilibrium/miller_batch.hpp"
#include "gfile_fixture.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
namespace fs = std::filesystem;

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double a, double b, double tolerance) {
    check(std::isfinite(a) && std::abs(a - b) <= tolerance,
          "unexpected Miller shape or scan coordinate");
}
std::string rejected(const auto& function) {
    try {
        function();
    } catch (const std::exception& e) { return e.what(); }
    throw std::runtime_error("invalid scan was accepted");
}
struct Row {
    std::string file;
    std::array<double, 13> values;
};
std::vector<Row> read_index(const fs::path& path) {
    std::ifstream input(path / "index.csv");
    std::string line;
    check(bool(std::getline(input, line)), "missing index");
    std::vector<Row> result;
    while (std::getline(input, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        auto& entry = result.emplace_back();
        row >> entry.file;
        for (auto& x : entry.values) check(bool(row >> x), "invalid index row");
        check(fs::is_regular_file(path / entry.file), "missing input file");
    }
    return result;
}

void exercise(const fs::path& work) {
    const auto source = work / "source.g";
    gfile_test::Settings fixture;
    fixture.axis_z = .25;
    fixture.field_sign = -1;
    std::ofstream(source) << gfile_test::fixture(fixture);
    using equilibrium::MillerRange;
    equilibrium::MillerScan scan;
    scan.r0 = MillerRange{6.0, 6.1, .1};
    scan.kappa_upper = MillerRange{1.4, 1.65, .2};
    scan.kappa_lower = MillerRange{1.4, 1.2, -.2};
    scan.delta_upper = MillerRange{0, .2, .2};
    scan.delta_lower = MillerRange{-.2, 0, .2};
    equilibrium::GfileInputSettings settings;
    settings.mpol = 16;
    const auto batch = work / "scan";
    check(equilibrium::write_miller_batch(source, batch, scan, settings) == 32,
          "wrong Cartesian count");
    auto rows = read_index(batch);
    check(rows.size() == 32, "wrong index length");
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& x = rows[i].values;
        near(x[0], 6.0 + .1 * (i / 16), 1e-14);
        near(x[1], 1.4 + .2 * ((i / 8) % 2), 1e-14);
        near(x[2], 1.4 - .2 * ((i / 4) % 2), 1e-14);
        near(x[3], .2 * ((i / 2) % 2), 1e-14);
        near(x[4], -.2 + .2 * (i % 2), 1e-14);
        near(x[5], gfile_test::A, 1e-8);
        near(x[6], .25, 1e-12);
        near(x[7], x[0], 1e-3 * x[5]);
        near(x[8], x[5], 1e-3 * x[5]);
        for (int j = 0; j < 4; ++j) near(x[j + 9], x[j + 1], 1e-3);
    }
    const auto original = equilibrium::gfile_input_to_json(
        equilibrium::gfile_to_cumes_input(source, settings));
    std::istringstream template_lines(original);
    std::ifstream generated(batch / rows.back().file);
    const std::string json((std::istreambuf_iterator<char>(generated)), {});
    std::string line;
    while (std::getline(template_lines, line))
        if (line.find("\"am\"") != line.npos ||
            line.find("\"ai\"") != line.npos ||
            line.find("\"phiedge\"") != line.npos ||
            line.find("\"niter_array\"") != line.npos)
            check(json.find(line) != json.npos,
                  "source pressure/iota/flux/schedule changed");
    rejected([&] {
        equilibrium::write_miller_batch(source, batch, scan, settings);
    });
    check(read_index(batch).size() == 32, "overwrote completed batch");

    // Omitted ranges use the source polygon even if a source Fourier fit would
    // fail.
    fixture.nonstar = true;
    std::ofstream(source) << gfile_test::fixture(fixture);
    const auto base =
        equilibrium::calculate_lcfs_shape(equilibrium::read_gfile_lcfs(source));
    const auto one = work / "defaults";
    check(equilibrium::write_miller_batch(source, one, {}, settings) == 1,
          "omitted ranges failed");
    const auto x = read_index(one).at(0).values;
    near(x[0], base.r_center, 1e-12);
    near(x[1], base.elongation_upper, 1e-12);
    near(x[2], base.elongation_lower, 1e-12);
    near(x[3], base.triangularity_upper, 1e-12);
    near(x[4], base.triangularity_lower, 1e-12);

    const auto bad = work / "failed";
    const auto fail = [&] {
        return rejected([&] {
            equilibrium::write_miller_batch(source, bad, scan, settings);
        });
    };
    scan.max_cases = 31;
    fail();
    check(!fs::exists(bad), "oversized product created output");
    scan = {};
    scan.r0 = MillerRange{1, 2, 0};
    fail();
    scan.r0 = MillerRange{2, 1, .1};
    fail();
    scan.r0 = MillerRange{1, 2, std::numeric_limits<double>::infinity()};
    fail();
    scan.r0 = MillerRange{6, 6.00000000000001, 1e-16};
    fail();
    scan.r0 = MillerRange{.5, .5, 1};
    fail();
    scan = {};
    scan.delta_lower = MillerRange{1, 1, 1};
    fail();
    scan = {};
    scan.kappa_upper = MillerRange{0, 0, 1};
    fail();
    // At mpol=2 an ellipse is exact. The second case needs additional modes.
    scan = {};
    scan.kappa_upper = scan.kappa_lower = MillerRange{1.4, 1.4, 1};
    scan.delta_upper = MillerRange{0, .5, .5};
    scan.delta_lower = MillerRange{0, 0, 1};
    settings.mpol = 2;
    check(fail().find("Miller case 1 ") != std::string::npos,
          "wrong failed scan case");
    check(!fs::exists(bad), "failed fit left a partial batch");
}
}  // namespace
int main(int argc, char** argv) {
    try {
        if (argc != 2)
            throw std::invalid_argument("expected fixture directory");
        // Use a new directory on every run without touching other test files.
        auto work = std::filesystem::path(argv[1]) / "miller-regression";
        if (!std::filesystem::create_directory(work))
            throw std::runtime_error("test directory already exists");
        try {
            exercise(work);
        } catch (...) {
            std::filesystem::remove_all(work);
            throw;
        }
        std::filesystem::remove_all(work);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
