#include "equilibrium/boozer_binary.hpp"
#include "equilibrium/boozer_output.hpp"
#include "equilibrium/io.hpp"
#include "hdf5_failure.hpp"

#if defined(__linux__) && !defined(__EMSCRIPTEN__)
#include <csignal>

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#ifdef EQUILIBRIUM_HAVE_NETCDF
#include <netcdf.h>
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
#include <hdf5.h>
#endif

#include <cmath>
#include <complex>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void expect(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

equilibrium::BoozerResult make_result() {
    equilibrium::BoozerResult original;
    original.source_format_version = 8;
    original.source_ns = 3;
    original.source_ntheta = 8;
    original.source_nzeta = 6;
    original.source_mpol = 3;
    original.source_ntor = 2;
    original.nfp = 5;
    original.radial_order = equilibrium::RadialInterpolationOrder::TWO_POINT;
    original.resonance_tolerance = 2.0e-11;
    original.s = {0.5, 1.0};
    original.iota = {0.4, 0.5};
    original.grid.source_ns = 3;
    original.grid.first_surface = 1;
    original.grid.ntheta = 2;
    original.grid.nzeta = 2;
    original.grid.b = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    original.grid.sqrtg_b = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0};
    original.grid.b2j00 = {-10.0, -11.0};
    original.spectrum.source_ns = 3;
    original.spectrum.first_surface = 1;
    original.spectrum.mmax = 0;
    original.spectrum.nmax = 0;
    original.spectrum.m.resize(1, 0);
    original.spectrum.n.resize(1, 0);
    original.spectrum.r = {{1.0, 0.0}, {3.0, 0.0}};
    original.spectrum.z = {{0.0, 0.0}, {0.0, 0.0}};
    original.spectrum.nu = {{0.0, 0.0}, {0.0, 0.0}};
    return original;
}

equilibrium::BoozerResult make_asymmetric_result() {
    auto result = make_result();
    result.source_format_version = 10;
    result.lasym = true;
    result.grid.ntheta = 16;
    result.grid.nzeta = 12;
    constexpr std::size_t POINTS = 16 * 12;
    result.grid.r.resize(2 * POINTS);
    result.grid.z.resize(2 * POINTS);
    result.grid.nu.resize(2 * POINTS);
    result.grid.b.assign(2 * POINTS, 2.0);
    result.grid.sqrtg_b.assign(2 * POINTS, -1.0);
    for (std::size_t surface = 0; surface < 2; ++surface) {
        const double scale = static_cast<double>(surface + 1);
        for (int n = 0; n < 12; ++n) {
            const double alpha = 2.0 * std::numbers::pi * n / 12.0;
            for (int m = 0; m < 16; ++m) {
                const double theta = 2.0 * std::numbers::pi * m / 16.0;
                const double value =
                    scale *
                    (3.0 + 0.2 * std::sin(theta) + 0.3 * std::sin(alpha) +
                     0.4 * std::sin(theta) * std::cos(alpha) +
                     0.5 * std::cos(theta) * std::sin(alpha) +
                     0.6 * std::sin(theta) * std::sin(alpha) +
                     0.7 * std::cos(theta) * std::cos(alpha));
                const auto index =
                    surface * POINTS + static_cast<std::size_t>(n * 16 + m);
                result.grid.r[index] = value;
                result.grid.z[index] = -2.0 * value;
                result.grid.nu[index] = 0.01 * value;
            }
        }
    }
    auto& spectrum = result.spectrum;
    spectrum.mmax = spectrum.nmax = 1;
    spectrum.m.clear();
    spectrum.n.clear();
    spectrum.r.clear();
    spectrum.z.clear();
    spectrum.nu.clear();
    for (int n = -1; n <= 1; ++n)
        for (int m = -1; m <= 1; ++m) {
            spectrum.m.push_back(m);
            spectrum.n.push_back(n);
        }
    // Analytic integral coefficients of the manufactured Fourier polynomial.
    for (int surface = 0; surface < 2; ++surface)
        for (int n = -1; n <= 1; ++n)
            for (int m = -1; m <= 1; ++m) {
                std::complex<double> value;
                if (m == 0 && n == 0)
                    value = 3.0;
                else if (m == 0)
                    value = {0.0, -0.15 * n};
                else if (n == 0)
                    value = {0.0, -0.1 * m};
                else
                    value = {(0.7 - 0.6 * m * n) / 4.0,
                             -(0.4 * m + 0.5 * n) / 4.0};
                value *=
                    (surface + 1) * 4.0 * std::numbers::pi * std::numbers::pi;
                spectrum.r.push_back(value);
                spectrum.z.push_back(-2.0 * value);
                spectrum.nu.push_back(0.01 * value);
            }
    return result;
}

equilibrium::BoozerResult with_tracing(equilibrium::BoozerResult result) {
    auto& tracing = result.tracing.emplace();
    tracing.phiedge = 7.5;
    tracing.aphi = {1.0, 0.2};
    tracing.toroidal_flux_derivative = {-2.0, -3.0, -4.0};
    tracing.poloidal_flux_derivative = {-0.6, -1.2, -2.0};
    tracing.axis_r_cos = {3.0, 0.1, 0.05};
    tracing.axis_r_sin = {0.0, result.lasym ? 0.02 : 0.0, 0.0};
    tracing.axis_z_cos = {result.lasym ? 0.12 : 0.0, 0.0, 0.0};
    tracing.axis_z_sin = {0.0, 0.03, 0.01};
    tracing.axis_half_s = {0.25, 0.75};
    for (std::size_t field = 0; field < tracing.axis_half_fields.size();
         ++field) {
        auto& values = tracing.axis_half_fields[field];
        values.resize(2 * 8 * 6);
        for (std::size_t i = 0; i < values.size(); ++i)
            values[i] = -2.0 - 0.1 * static_cast<double>(field) -
                        0.001 * static_cast<double>(i);
    }
    return result;
}

void with_magnetic(equilibrium::BoozerResult& result) {
    auto& magnetic = result.tracing->magnetic.emplace();
    for (std::size_t field = 0; field < magnetic.fields.size(); ++field)
        magnetic.fields[field].assign(result.grid.b.size(),
                                      -0.3 * static_cast<double>(field + 1));
    result.grid.source_theta_shift.resize(result.grid.b.size());
    result.grid.nu.resize(result.grid.b.size());
    for (std::size_t i = 0; i < result.grid.b.size(); ++i) {
        // Deliberately beyond the truncated spectrum and outside [-pi,pi].
        // The component map must preserve these samples without wrapping.
        result.grid.source_theta_shift[i] =
            4.0 + 0.03 * std::sin(static_cast<double>(i));
        result.grid.nu[i] = 0.07 * std::cos(0.3 * static_cast<double>(i));
    }
}

void with_profiles(equilibrium::BoozerResult& result) {
    auto& profiles = result.profiles.emplace();
    profiles.pressure = {1200.0, 100.0};
    profiles.current_i = {0.2, 0.3};
    profiles.current_g = {-1.7, -2.2};
}

using NamedArray = std::pair<const char*, std::span<const double>>;
std::vector<NamedArray> tracing_arrays(
    const equilibrium::BoozerResult& result) {
    const auto& data = *result.tracing;
    std::vector<NamedArray> arrays = {
        {"s", result.s},
        {"iota", result.iota},
        {"aphi", data.aphi},
        {"toroidal_flux_derivative", data.toroidal_flux_derivative},
        {"poloidal_flux_derivative", data.poloidal_flux_derivative},
        {"axis_r_cos", data.axis_r_cos},
        {"axis_r_sin", data.axis_r_sin},
        {"axis_z_cos", data.axis_z_cos},
        {"axis_z_sin", data.axis_z_sin},
        {"axis_half_s", data.axis_half_s}};
    const char* names[] = {"axis_half_sqrtg", "axis_half_bsups",
                           "axis_half_bsupu", "axis_half_bsupv",
                           "axis_half_bsubs", "axis_half_bsubu",
                           "axis_half_bsubv"};
    for (std::size_t i = 0; i < data.axis_half_fields.size(); ++i)
        arrays.emplace_back(names[i], data.axis_half_fields[i]);
    if (data.magnetic) {
        const auto& magnetic = *data.magnetic;
        const char* names[] = {"boozer_sqrtg", "boozer_bsups", "boozer_bsupu",
                               "boozer_bsupv", "boozer_bsubs", "boozer_bsubu",
                               "boozer_bsubv"};
        for (std::size_t i = 0; i < magnetic.fields.size(); ++i)
            arrays.emplace_back(names[i], magnetic.fields[i]);
        arrays.emplace_back("source_theta_shift",
                            result.grid.source_theta_shift);
        arrays.emplace_back("nu_grid", result.grid.nu);
    }
    if (result.profiles) {
        arrays.emplace_back("pressure", result.profiles->pressure);
        arrays.emplace_back("current_i", result.profiles->current_i);
        arrays.emplace_back("current_g", result.profiles->current_g);
    }
    return arrays;
}

void test_tracing_round_trip(bool lasym,
                             bool magnetic_fields,
                             bool surface_profiles = false) {
    auto original =
        with_tracing(lasym ? make_asymmetric_result() : make_result());
    if (magnetic_fields) with_magnetic(original);
    if (surface_profiles) with_profiles(original);
    const auto path = std::filesystem::current_path() /
                      "magnetic-coordinate-boozer-tracing.bin";
    auto legacy = original;
    legacy.tracing.reset();
    legacy.profiles.reset();
    equilibrium::write_boozer_binary(path, legacy, "manufactured");
    const auto reference = equilibrium::read_boozer_binary(path);
    equilibrium::write_boozer_binary(path, original, "manufactured");
    {
        std::ifstream input(path, std::ios::binary);
        std::string magic(8, '\0');
        input.read(magic.data(), 8);
        expect(magic == (surface_profiles  ? "MCBOOZ07"
                         : magnetic_fields ? "MCBOOZ06"
                                           : "MCBOOZ05"),
               "incorrect tracing format version");
    }
    const auto loaded = equilibrium::read_boozer_binary(path);
    expect(loaded.result.lasym == lasym && loaded.result.tracing.has_value(),
           "tracing metadata");
    expect(loaded.result.tracing->magnetic.has_value() == magnetic_fields,
           "legacy v5 must not invent Boozer components or maps");
    expect(loaded.result.profiles.has_value() == surface_profiles &&
               !reference.result.profiles,
           "legacy files must not invent surface profiles");
    const auto expected = tracing_arrays(original);
    const auto actual = tracing_arrays(loaded.result);
    expect(loaded.result.tracing->phiedge == original.tracing->phiedge,
           "signed input flux");
    for (std::size_t i = 0; i < expected.size(); ++i)
        expect(std::vector<double>(actual[i].second.begin(),
                                   actual[i].second.end()) ==
                   std::vector<double>(expected[i].second.begin(),
                                       expected[i].second.end()),
               expected[i].first);
    expect(loaded.result.grid.b == original.grid.b &&
               loaded.result.grid.sqrtg_b == original.grid.sqrtg_b,
           "v5 must preserve the non-axis grids");
    expect(loaded.result.spectrum.r == reference.result.spectrum.r &&
               loaded.result.spectrum.z == reference.result.spectrum.z &&
               loaded.result.spectrum.nu == reference.result.spectrum.nu,
           "v5 must preserve spectra");
    std::filesystem::resize_file(path, std::filesystem::file_size(path) - 1);
    bool truncated = false;
    try {
        static_cast<void>(equilibrium::read_boozer_binary(path));
    } catch (const std::runtime_error&) { truncated = true; }
    expect(truncated, "truncated tracing data accepted");
    if (surface_profiles) {
        equilibrium::write_boozer_binary(path, original, "manufactured");
        {
            std::fstream output(
                path, std::ios::binary | std::ios::in | std::ios::out);
            output.seekp(-8, std::ios::end);
            const unsigned char nan[] = {0, 0, 0, 0, 0, 0, 0xf8, 0x7f};
            output.write(reinterpret_cast<const char*>(nan), 8);
        }
        bool rejected = false;
        try {
            static_cast<void>(equilibrium::read_boozer_binary(path));
        } catch (const std::invalid_argument&) { rejected = true; }
        expect(rejected, "non-finite v7 profile accepted by reader");
    }
    for (int invalid_case = 0; invalid_case < (surface_profiles  ? 20
                                               : magnetic_fields ? 11
                                                                 : 6);
         ++invalid_case) {
        auto invalid = original;
        auto& data = *invalid.tracing;
        switch (invalid_case) {
            case 0:
                data.toroidal_flux_derivative[0] = 0.0;
                break;
            case 1:
                data.poloidal_flux_derivative[1] *= 2.0;
                break;
            case 2:
                data.axis_half_s[1] = data.axis_half_s[0];
                break;
            case 3:
                data.axis_r_cos.pop_back();
                break;
            case 4:
                data.axis_half_fields[6][0] =
                    std::numeric_limits<double>::quiet_NaN();
                break;
            case 5:
                data.aphi.assign(1, 0.0);
                break;
            case 6:
                data.magnetic->fields[0][0] = 0.0;
                break;
            case 7:
                data.magnetic->fields[3].pop_back();
                break;
            case 8:
                data.magnetic->fields[1][0] =
                    std::numeric_limits<double>::quiet_NaN();
                break;
            case 9:
                invalid.grid.source_theta_shift[0] =
                    std::numeric_limits<double>::infinity();
                break;
            case 10:
                invalid.grid.nu.pop_back();
                break;
            case 11:
                invalid.profiles->pressure.pop_back();
                break;
            case 12:
                invalid.profiles->current_i.pop_back();
                break;
            case 13:
                invalid.profiles->current_g.pop_back();
                break;
            case 14:
                invalid.profiles->pressure[0] =
                    std::numeric_limits<double>::quiet_NaN();
                break;
            case 15:
                invalid.profiles->current_i[0] =
                    std::numeric_limits<double>::infinity();
                break;
            case 16:
                invalid.profiles->current_g[0] =
                    std::numeric_limits<double>::quiet_NaN();
                break;
            case 17:
                invalid.profiles->current_g[0] =
                    -invalid.iota[0] * invalid.profiles->current_i[0];
                break;
            case 18:
                data.magnetic.reset();
                break;
            case 19:
                invalid.tracing.reset();
                break;
        }
        bool rejected = false;
        try {
            equilibrium::write_boozer_binary(path, invalid, "invalid");
        } catch (const std::invalid_argument&) { rejected = true; }
        expect(rejected, "invalid tracing data accepted");
    }
    std::filesystem::remove(path);
}

void test_flux_list_round_trip() {
    for (const bool lasym : {false, true}) {
        auto original =
            with_tracing(lasym ? make_asymmetric_result() : make_result());
        with_magnetic(original);
        with_profiles(original);
        original.source_format_version = 11;
        original.tracing->pflux_type = "list";
        original.tracing->aphi = {1.0, 2.0, 1.0};
        const auto path = std::filesystem::current_path() /
                          "magnetic-coordinate-flux-list.bin";
        equilibrium::write_boozer_binary(path, original, "flux-list");
        std::ifstream bytes(path, std::ios::binary);
        std::string magic(8, '\0');
        bytes.read(magic.data(), 8);
        bytes.close();
        expect(magic == "MCBOOZ08", "flux list requires Boozer version 8");
        const auto loaded = equilibrium::read_boozer_binary(path);
        expect(
            loaded.result.lasym == lasym &&
                loaded.result.source_format_version == 11 &&
                loaded.result.tracing->pflux_type == "list" &&
                loaded.result.tracing->aphi == original.tracing->aphi &&
                loaded.result.profiles->pressure == original.profiles->pressure,
            "flux list metadata and pressure round trip");
        std::filesystem::resize_file(path,
                                     std::filesystem::file_size(path) - 1);
        bool rejected = false;
        try {
            static_cast<void>(equilibrium::read_boozer_binary(path));
        } catch (const std::runtime_error&) { rejected = true; }
        expect(rejected, "truncated flux type rejected");
        original.tracing->aphi = {1, 0};
        rejected = false;
        try {
            equilibrium::write_boozer_binary(path, original, "invalid");
        } catch (const std::invalid_argument&) { rejected = true; }
        expect(rejected, "nonpositive flux derivative rejected");
        std::filesystem::remove(path);
    }
}

void test_asymmetric_round_trip() {
    auto original = make_asymmetric_result();
    const auto path = std::filesystem::current_path() /
                      "magnetic-coordinate-boozer-asymmetric.bin";
    equilibrium::write_boozer_binary(path, original, "asymmetric.bin");
    {
        std::ifstream input(path, std::ios::binary);
        std::string magic(8, '\0');
        input.read(magic.data(), 8);
        expect(magic == "MCBOOZ04", "asymmetric output requires version 4");
    }
    const auto loaded = equilibrium::read_boozer_binary(path);
    expect(loaded.result.lasym && loaded.result.source_format_version == 10,
           "asymmetric output metadata");
    expect(!loaded.result.tracing,
           "legacy v4 must not invent polarity or axis data");
    expect(!loaded.result.profiles, "legacy v4 must not invent profiles");
    const auto& spectrum = loaded.result.spectrum;
    const std::vector<double>* fields[] = {&original.grid.r, &original.grid.z,
                                           &original.grid.nu};
    const std::vector<std::complex<double>>* coefficients[] = {
        &spectrum.r, &spectrum.z, &spectrum.nu};
    for (int field = 0; field < 3; ++field) {
        for (std::size_t surface = 0; surface < 2; ++surface) {
            for (int n = 0; n < 12; ++n) {
                const double alpha = 2.0 * std::numbers::pi * n / 12.0;
                for (int m = 0; m < 16; ++m) {
                    const double theta = 2.0 * std::numbers::pi * m / 16.0;
                    std::complex<double> value{};
                    for (std::size_t mode = 0; mode < spectrum.m.size();
                         ++mode) {
                        value +=
                            (*coefficients[field])[surface * spectrum.m.size() +
                                                   mode] *
                            std::polar(1.0, spectrum.m[mode] * theta +
                                                spectrum.n[mode] * alpha);
                    }
                    value /= 4.0 * std::numbers::pi * std::numbers::pi;
                    const auto index = surface * 16 * 12 +
                                       static_cast<std::size_t>(n * 16 + m);
                    expect(std::abs(value - (*fields[field])[index]) < 3.0e-13,
                           "asymmetric real families must reconstruct the "
                           "analytic field");
                }
            }
        }
    }
    std::filesystem::resize_file(path, std::filesystem::file_size(path) - 1);
    bool truncated = false;
    try {
        static_cast<void>(equilibrium::read_boozer_binary(path));
    } catch (const std::runtime_error&) { truncated = true; }
    expect(truncated, "truncated asymmetric output accepted");
    const auto reject = [&] {
        try {
            equilibrium::write_boozer_binary(path, original, "invalid.bin");
        } catch (const std::invalid_argument&) { return; }
        throw std::runtime_error("invalid parity/Hermitian spectrum accepted");
    };
    original.lasym = false;
    reject();
    original.lasym = true;
    original.spectrum.r[0] += std::complex<double>(0.1, 0.2);
    reject();
    std::filesystem::remove(path);
}

void test_round_trip() {
    const auto original = make_result();

    const auto path = std::filesystem::current_path() /
                      "magnetic-coordinate-boozer-roundtrip.bin";
    equilibrium::write_boozer_binary(path, original, "source-output.bin");
    const auto loaded = equilibrium::read_boozer_binary(path);
    std::filesystem::remove(path);
    expect(loaded.source_path == "source-output.bin", "source provenance");
    expect(!loaded.result.tracing,
           "legacy v3 must not invent polarity or axis data");
    expect(!loaded.result.profiles, "legacy v3 must not invent profiles");
    expect(loaded.result.source_format_version == 8, "source version");
    expect(loaded.result.nfp == 5, "nfp");
    expect(loaded.result.resonance_tolerance == 2.0e-11, "resonance tolerance");
    expect(loaded.result.s == original.s, "surface coordinates");
    expect(loaded.result.iota == original.iota, "iota");
    expect(loaded.result.spectrum.r == original.spectrum.r, "R spectrum");
    expect(loaded.result.spectrum.z == original.spectrum.z, "Z spectrum");
    expect(loaded.result.spectrum.nu == original.spectrum.nu, "nu spectrum");
    expect(loaded.result.grid.b == original.grid.b, "B grid");
    expect(loaded.result.grid.sqrtg_b == original.grid.sqrtg_b,
           "Jacobian grid");
    expect(loaded.result.grid.b2j00 == original.grid.b2j00, "b2j00");
}

void test_nontrivial_real_parity_round_trip() {
    auto original = make_result();
    original.spectrum.mmax = 1;
    original.spectrum.nmax = 1;
    original.spectrum.m.clear();
    original.spectrum.n.clear();
    for (int n = -1; n <= 1; ++n) {
        for (int m = -1; m <= 1; ++m) {
            original.spectrum.m.push_back(m);
            original.spectrum.n.push_back(n);
        }
    }
    original.spectrum.r.assign(18, {});
    original.spectrum.z.assign(18, {});
    original.spectrum.nu.assign(18, {});
    const auto mode = [](int m, int n) {
        return static_cast<std::size_t>((n + 1) * 3 + (m + 1));
    };
    for (int surface = 0; surface < 2; ++surface) {
        const double scale = std::numbers::pi * std::numbers::pi *
                             static_cast<double>(surface + 1);
        const auto coefficient = [surface](std::size_t mode_index) {
            return static_cast<std::size_t>(surface * 9) + mode_index;
        };
        original.spectrum.r[coefficient(mode(1, 1))] = {1.5 * scale, 0.0};
        original.spectrum.r[coefficient(mode(1, -1))] = {2.5 * scale, 0.0};
        original.spectrum.r[coefficient(mode(-1, -1))] = {1.5 * scale, 0.0};
        original.spectrum.r[coefficient(mode(-1, 1))] = {2.5 * scale, 0.0};
        original.spectrum.z[coefficient(mode(1, 1))] = {0.0, -0.5 * scale};
        original.spectrum.z[coefficient(mode(1, -1))] = {0.0, -2.0 * scale};
        original.spectrum.z[coefficient(mode(-1, -1))] = {0.0, 0.5 * scale};
        original.spectrum.z[coefficient(mode(-1, 1))] = {0.0, 2.0 * scale};
        original.spectrum.nu[coefficient(mode(1, 1))] = {0.0, -0.5 * scale};
        original.spectrum.nu[coefficient(mode(1, -1))] = {0.0, 1.3 * scale};
        original.spectrum.nu[coefficient(mode(-1, -1))] = {0.0, 0.5 * scale};
        original.spectrum.nu[coefficient(mode(-1, 1))] = {0.0, -1.3 * scale};
    }

    const auto path = std::filesystem::current_path() /
                      "magnetic-coordinate-boozer-parity.bin";
    equilibrium::write_boozer_binary(path, original, "source.bin");
    const auto loaded = equilibrium::read_boozer_binary(path);
    std::filesystem::remove(path);
    const auto close = [](const auto& left, const auto& right) {
        if (left.size() != right.size()) return false;
        for (std::size_t index = 0; index < left.size(); ++index) {
            if (std::abs(left[index] - right[index]) > 1.0e-13) return false;
        }
        return true;
    };
    expect(close(loaded.result.spectrum.r, original.spectrum.r),
           "nontrivial R parity round trip");
    expect(close(loaded.result.spectrum.z, original.spectrum.z),
           "nontrivial Z parity round trip");
    expect(close(loaded.result.spectrum.nu, original.spectrum.nu),
           "nontrivial nu parity round trip");
}

void test_dispatch() {
    using equilibrium::BoozerOutputFormat;
    expect(equilibrium::resolve_boozer_output_spec("result.BIN").format ==
               BoozerOutputFormat::BINARY,
           "binary suffix dispatch");
    expect(equilibrium::resolve_boozer_output_spec("result.NC").format ==
               BoozerOutputFormat::NETCDF,
           "NetCDF suffix dispatch");
    expect(equilibrium::resolve_boozer_output_spec("result.HDF5").format ==
               BoozerOutputFormat::HDF5,
           "HDF5 suffix dispatch");
#ifndef EQUILIBRIUM_HAVE_NETCDF
    expect(!equilibrium::boozer_output_format_available(
               BoozerOutputFormat::NETCDF),
           "NetCDF unavailable without backend");
#endif
#ifndef EQUILIBRIUM_HAVE_HDF5
    expect(
        !equilibrium::boozer_output_format_available(BoozerOutputFormat::HDF5),
        "HDF5 unavailable without backend");
#endif
    try {
        static_cast<void>(
            equilibrium::resolve_boozer_output_spec("result.dat"));
    } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("unknown output suffix was accepted");
}

#ifdef EQUILIBRIUM_HAVE_NETCDF
void test_netcdf_output(bool lasym,
                        bool tracing,
                        bool magnetic_fields = false,
                        bool surface_profiles = false,
                        bool flux_list = false) {
    const auto path = std::filesystem::current_path() /
                      "magnetic-coordinate-boozer-output.nc";
    auto original = lasym ? make_asymmetric_result() : make_result();
    if (tracing) original = with_tracing(std::move(original));
    if (magnetic_fields) with_magnetic(original);
    if (surface_profiles) with_profiles(original);
    if (flux_list) {
        original.tracing->pflux_type = "list";
        original.tracing->aphi = {1, 2, 1};
    }
    equilibrium::write_boozer_output(path, original, "source-output.bin");
    int file = -1;
    expect(nc_open(path.string().c_str(), NC_NOWRITE, &file) == NC_NOERR,
           "open NetCDF output");
    int variable = -1;
    expect(nc_inq_varid(file, "rmncc", &variable) == NC_NOERR,
           "NetCDF rmncc exists");
    nc_type type = NC_NAT;
    int dimensions = 0;
    expect(nc_inq_var(file, variable, nullptr, &type, &dimensions, nullptr,
                      nullptr) == NC_NOERR &&
               type == NC_DOUBLE && dimensions == 2,
           "NetCDF rmncc is a real rank-2 array");
    expect(nc_inq_varid(file, "R_spectrum", &variable) == NC_ENOTVAR,
           "NetCDF has no complex spectrum variable");
    int alpha_dimension = -1;
    expect(nc_inq_dimid(file, "alpha", &alpha_dimension) == NC_NOERR,
           "NetCDF has field-period alpha dimension");
    int obsolete_zeta_dimension = -1;
    expect(nc_inq_dimid(file, "zeta", &obsolete_zeta_dimension) == NC_EBADDIM,
           "NetCDF has no ambiguous zeta dimension");
    int complex_dimension = -1;
    expect(nc_inq_dimid(file, "complex_component", &complex_dimension) ==
               NC_EBADDIM,
           "NetCDF has no complex dimension");
    expect(nc_inq_varid(file, "rmncc", &variable) == NC_NOERR,
           "find NetCDF rmncc");
    std::vector<double> rmncc(lasym ? 8 : 2);
    expect(nc_get_var_double(file, variable, rmncc.data()) == NC_NOERR,
           "read NetCDF rmncc");
    expect(std::abs(rmncc[0] - (lasym ? 3.0
                                      : 1.0 / (4.0 * std::numbers::pi *
                                               std::numbers::pi))) < 1.0e-15,
           "NetCDF real coefficient normalization");
    if (lasym) {
        int flag = 0;
        expect(nc_get_att_int(file, NC_GLOBAL, "lasym", &flag) == NC_NOERR &&
                   flag == 1,
               "NetCDF asymmetric metadata");
        const char* names[] = {"rmnsc", "rmncs",  "zmncc",
                               "zmnss", "numncc", "numnss"};
        const double expected[] = {0.4, 0.5, -1.4, -1.2, 0.007, 0.006};
        for (int family = 0; family < 6; ++family) {
            expect(nc_inq_varid(file, names[family], &variable) == NC_NOERR,
                   "NetCDF complementary family missing");
            std::vector<double> values(8);
            expect(nc_get_var_double(file, variable, values.data()) == NC_NOERR,
                   "read NetCDF complementary family");
            expect(std::abs(values[3] - expected[family]) < 1.0e-14 &&
                       std::abs(values[7] - 2.0 * expected[family]) < 2.0e-14,
                   "NetCDF complementary family normalization");
        }
    }
    if (tracing) {
        double phiedge = 0.0;
        expect(nc_get_att_double(file, NC_GLOBAL, "phiedge", &phiedge) ==
                       NC_NOERR &&
                   phiedge == 7.5,
               "NetCDF signed flux");
        char schema[64] = {};
        expect(nc_get_att_text(file, NC_GLOBAL, "schema", schema) == NC_NOERR &&
                   std::string(schema) ==
                       (flux_list          ? "magnetic-coordinate-boozer-v8"
                        : surface_profiles ? "magnetic-coordinate-boozer-v7"
                        : magnetic_fields  ? "magnetic-coordinate-boozer-v6"
                                           : "magnetic-coordinate-boozer-v5"),
               "NetCDF schema");
        if (flux_list) {
            char type[16] = {};
            expect(nc_get_att_text(file, NC_GLOBAL, "pflux_type", type) ==
                           NC_NOERR &&
                       std::string(type) == "list",
                   "NetCDF flux selector");
        }
        for (const auto& [name, expected] : tracing_arrays(original)) {
            expect(nc_inq_varid(file, name, &variable) == NC_NOERR, name);
            std::vector<double> actual(expected.size());
            expect(
                nc_get_var_double(file, variable, actual.data()) == NC_NOERR &&
                    actual ==
                        std::vector<double>(expected.begin(), expected.end()),
                name);
            if (surface_profiles && (std::string(name) == "pressure" ||
                                     std::string(name) == "current_i" ||
                                     std::string(name) == "current_g")) {
                char units[8] = {};
                expect(nc_get_att_text(file, variable, "units", units) ==
                               NC_NOERR &&
                           std::string(units) ==
                               (std::string(name) == "pressure" ? "Pa" : "T m"),
                       "NetCDF profile units");
                int rank = 0, dim = -1;
                std::size_t size = 0;
                expect(nc_inq_varndims(file, variable, &rank) == NC_NOERR &&
                           rank == 1 &&
                           nc_inq_vardimid(file, variable, &dim) == NC_NOERR &&
                           nc_inq_dimlen(file, dim, &size) == NC_NOERR &&
                           size == original.s.size(),
                       "NetCDF profile surface grid");
            }
            if (std::string(name) == "axis_half_bsupv") {
                int dims[3] = {};
                expect(nc_inq_vardimid(file, variable, dims) == NC_NOERR,
                       "NetCDF native axis grid");
                for (int i = 0; i < 3; ++i) {
                    std::size_t count = 0;
                    expect(nc_inq_dimlen(file, dims[i], &count) == NC_NOERR &&
                               count == static_cast<std::size_t>(i == 0   ? 2
                                                                 : i == 1 ? 6
                                                                          : 8),
                           "NetCDF axis grid dimensions");
                }
            }
        }
    }
    if (!surface_profiles)
        expect(nc_inq_varid(file, "pressure", &variable) == NC_ENOTVAR,
               "legacy NetCDF must not invent pressure");
    expect(nc_close(file) == NC_NOERR, "close NetCDF output");
    std::filesystem::remove(path);
}
#endif

#ifdef EQUILIBRIUM_HAVE_HDF5
void test_hdf5_output(bool lasym,
                      bool tracing,
                      bool magnetic_fields = false,
                      bool surface_profiles = false,
                      bool flux_list = false) {
    const auto path = std::filesystem::current_path() /
                      "magnetic-coordinate-boozer-output.h5";
    auto original = lasym ? make_asymmetric_result() : make_result();
    if (tracing) original = with_tracing(std::move(original));
    if (magnetic_fields) with_magnetic(original);
    if (surface_profiles) with_profiles(original);
    if (flux_list) {
        original.tracing->pflux_type = "list";
        original.tracing->aphi = {1, 2, 1};
    }
    equilibrium::write_boozer_output(path, original, "source-output.bin");
    const hid_t file =
        H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    expect(file >= 0, "open HDF5 output");
    if (surface_profiles) {
        const hid_t attribute = H5Aopen(file, "schema", H5P_DEFAULT);
        const hid_t schema_type = H5Aget_type(attribute);
        char schema[64] = {};
        expect(attribute >= 0 && schema_type >= 0 &&
                   H5Aread(attribute, schema_type, schema) >= 0 &&
                   std::string(schema) ==
                       (flux_list ? "magnetic-coordinate-boozer-v8"
                                  : "magnetic-coordinate-boozer-v7"),
               "HDF5 v7 schema");
        H5Tclose(schema_type);
        H5Aclose(attribute);
    }
    if (flux_list) {
        const hid_t attribute = H5Aopen(file, "pflux_type", H5P_DEFAULT);
        const hid_t type = H5Aget_type(attribute);
        char value[16] = {};
        expect(attribute >= 0 && type >= 0 &&
                   H5Aread(attribute, type, value) >= 0 &&
                   std::string(value) == "list",
               "HDF5 flux selector");
        H5Tclose(type);
        H5Aclose(attribute);
    }
    const hid_t dataset = H5Dopen2(file, "rmncc", H5P_DEFAULT);
    expect(dataset >= 0, "HDF5 rmncc exists");
    const hid_t space = H5Dget_space(dataset);
    expect(space >= 0 && H5Sget_simple_extent_ndims(space) == 2,
           "HDF5 rmncc is rank 2");
    const hid_t type = H5Dget_type(dataset);
    expect(type >= 0 && H5Tget_class(type) == H5T_FLOAT,
           "HDF5 rmncc is real-valued");
    expect(H5Lexists(file, "R_spectrum", H5P_DEFAULT) == 0,
           "HDF5 has no complex spectrum dataset");
    expect(H5Lexists(file, "alpha", H5P_DEFAULT) > 0,
           "HDF5 has field-period alpha dataset");
    expect(H5Lexists(file, "zeta", H5P_DEFAULT) == 0,
           "HDF5 has no ambiguous zeta dataset");
    double rmncc[8] = {};
    expect(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   rmncc) >= 0,
           "read HDF5 rmncc");
    expect(std::abs(rmncc[0] - (lasym ? 3.0
                                      : 1.0 / (4.0 * std::numbers::pi *
                                               std::numbers::pi))) < 1.0e-15,
           "HDF5 real coefficient normalization");
    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(dataset);
    if (lasym) {
        const hid_t attribute = H5Aopen(file, "lasym", H5P_DEFAULT);
        int flag = 0;
        expect(attribute >= 0 &&
                   H5Aread(attribute, H5T_NATIVE_INT, &flag) >= 0 && flag == 1,
               "HDF5 asymmetric metadata");
        H5Aclose(attribute);
        const char* names[] = {"rmnsc", "rmncs",  "zmncc",
                               "zmnss", "numncc", "numnss"};
        const double expected[] = {0.4, 0.5, -1.4, -1.2, 0.007, 0.006};
        for (int family = 0; family < 6; ++family) {
            const hid_t extra = H5Dopen2(file, names[family], H5P_DEFAULT);
            expect(extra >= 0, "HDF5 complementary family missing");
            double values[8] = {};
            expect(H5Dread(extra, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, values) >= 0,
                   "read HDF5 complementary family");
            expect(std::abs(values[3] - expected[family]) < 1.0e-14 &&
                       std::abs(values[7] - 2.0 * expected[family]) < 2.0e-14,
                   "HDF5 complementary family normalization");
            H5Dclose(extra);
        }
    }
    if (tracing) {
        const hid_t flux_attribute = H5Aopen(file, "phiedge", H5P_DEFAULT);
        double phiedge = 0.0;
        expect(flux_attribute >= 0 &&
                   H5Aread(flux_attribute, H5T_NATIVE_DOUBLE, &phiedge) >= 0 &&
                   phiedge == 7.5,
               "HDF5 signed flux");
        H5Aclose(flux_attribute);
        for (const auto& [name, expected] : tracing_arrays(original)) {
            const hid_t array = H5Dopen2(file, name, H5P_DEFAULT);
            expect(array >= 0, name);
            if (surface_profiles && (std::string(name) == "pressure" ||
                                     std::string(name) == "current_i" ||
                                     std::string(name) == "current_g")) {
                const hid_t attribute = H5Aopen(array, "units", H5P_DEFAULT);
                const hid_t units_type = H5Aget_type(attribute);
                char units[8] = {};
                expect(attribute >= 0 && units_type >= 0 &&
                           H5Aread(attribute, units_type, units) >= 0 &&
                           std::string(units) ==
                               (std::string(name) == "pressure" ? "Pa" : "T m"),
                       "HDF5 profile units");
                H5Tclose(units_type);
                H5Aclose(attribute);
                const hid_t profile_space = H5Dget_space(array);
                hsize_t size = 0;
                expect(profile_space >= 0 &&
                           H5Sget_simple_extent_ndims(profile_space) == 1 &&
                           H5Sget_simple_extent_dims(profile_space, &size,
                                                     nullptr) == 1 &&
                           size == original.s.size(),
                       "HDF5 profile surface grid");
                H5Sclose(profile_space);
            }
            std::vector<double> actual(expected.size());
            expect(H5Dread(array, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, actual.data()) >= 0 &&
                       actual == std::vector<double>(expected.begin(),
                                                     expected.end()),
                   name);
            if (std::string(name) == "axis_half_bsupv") {
                const hid_t native_space = H5Dget_space(array);
                hsize_t dims[3] = {};
                expect(H5Sget_simple_extent_dims(native_space, dims, nullptr) ==
                               3 &&
                           dims[0] == 2 && dims[1] == 6 && dims[2] == 8,
                       "HDF5 axis grid dimensions");
                H5Sclose(native_space);
            }
            H5Dclose(array);
        }
    }
    if (!surface_profiles)
        expect(H5Lexists(file, "pressure", H5P_DEFAULT) == 0,
               "legacy HDF5 must not invent pressure");
    H5Fclose(file);
    std::filesystem::remove(path);
}
#endif

// Failures before publication must preserve the previous file and remove only
// the writer's own temporary file, including a failed library create call.
void test_publication_and_conversion() {
    const auto result = make_result();
    const std::filesystem::path input = "boozer-convert-source.bin";
    equilibrium::write_boozer_binary(input, result, "original-source.bin");
    for (const auto format : {equilibrium::BoozerOutputFormat::BINARY,
                              equilibrium::BoozerOutputFormat::NETCDF,
                              equilibrium::BoozerOutputFormat::HDF5}) {
        if (!equilibrium::boozer_output_format_available(format)) continue;
        const std::string path =
            "boozer-publication" +
            std::string(equilibrium::boozer_output_suffix(format));
        const std::string sentinel = "previous destination";
        auto seed = [&](const std::string& name) {
            std::ofstream file(name, std::ios::binary);
            file << sentinel;
        };
        auto preserved = [&](const std::string& name) {
            std::ifstream file(name, std::ios::binary);
            return std::string(std::istreambuf_iterator<char>(file), {}) ==
                   sentinel;
        };
        auto no_temporary = [&] {
            for (const auto& entry : std::filesystem::directory_iterator("."))
                expect(!entry.path().filename().string().starts_with(path +
                                                                     ".tmp."),
                       "failed writer left a temporary file");
        };
        seed(path);
        auto invalid = result;
        invalid.s.clear();
        bool rejected = false;
        try {
            equilibrium::write_boozer_output(path, invalid, input);
        } catch (const std::exception&) { rejected = true; }
        expect(rejected && preserved(path),
               "invalid result replaced destination");
        no_temporary();
#if defined(__linux__) && !defined(__EMSCRIPTEN__)
#ifdef EQUILIBRIUM_HAVE_HDF5
        if (format == equilibrium::BoozerOutputFormat::HDF5) {
            equilibrium::test::hdf5_write_failures(
                path,
                [&](const std::string& destination) {
                    equilibrium::write_boozer_output(destination, result,
                                                     input);
                },
                [&](const std::string& destination) {
                    const hid_t file = H5Fopen(destination.c_str(),
                                               H5F_ACC_RDONLY, H5P_DEFAULT);
                    expect(file >= 0, "cannot reopen recovered Boozer file");
                    const hid_t dataset = H5Dopen2(file, "B", H5P_DEFAULT);
                    expect(dataset >= 0, "missing recovered Boozer field");
                    std::vector<double> values(result.grid.b.size());
                    const auto read =
                        H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, values.data());
                    H5Dclose(dataset);
                    H5Fclose(file);
                    expect(read >= 0 && values == result.grid.b,
                           "recovered Boozer field changed");
                });
        } else
#endif
        {
            const pid_t child = fork();
            expect(child >= 0, "cannot fork write-failure test");
            if (child == 0) {
                std::signal(SIGXFSZ, SIG_IGN);
                const rlimit limit{256, 256};
                if (setrlimit(RLIMIT_FSIZE, &limit) != 0) _exit(2);
                try {
                    equilibrium::write_boozer_output(path, result, input);
                } catch (const std::exception&) { _exit(0); }
                _exit(1);
            }
            int status = 0;
            expect(waitpid(child, &status, 0) == child && WIFEXITED(status) &&
                       WEXITSTATUS(status) == 0,
                   "write failure did not propagate");
            expect(preserved(path), "write failure replaced destination");
            no_temporary();
        }
#endif
        std::filesystem::remove(path);
        std::filesystem::create_directory(path);
        seed(path + "/keep");
        rejected = false;
        try {
            equilibrium::write_boozer_output(path, result, input);
        } catch (const std::exception&) { rejected = true; }
        expect(rejected && preserved(path + "/keep"),
               "rename failure damaged destination");
        no_temporary();
        std::filesystem::remove_all(path);
        seed(path + ".tmp");
        expect(equilibrium::convert(input.string(), path).has_value(),
               "Boozer conversion failed");
        expect(preserved(path + ".tmp"),
               "writer removed another temporary file");
        no_temporary();
        if (format == equilibrium::BoozerOutputFormat::BINARY) {
            const auto converted = equilibrium::read_boozer_binary(path);
            expect(converted.source_path == "original-source.bin" &&
                       converted.result.grid.b == result.grid.b,
                   "conversion lost source path or field values");
        }
        std::filesystem::remove(path);
        std::filesystem::remove(path + ".tmp");
    }
    std::filesystem::remove(input);
}

}  // namespace

int main() {
    try {
        test_round_trip();
        test_flux_list_round_trip();
        for (bool magnetic_fields : {false, true}) {
            test_tracing_round_trip(false, magnetic_fields);
            test_tracing_round_trip(true, magnetic_fields);
        }
        test_tracing_round_trip(false, true, true);
        test_tracing_round_trip(true, true, true);
        test_asymmetric_round_trip();
        test_nontrivial_real_parity_round_trip();
        test_dispatch();
        test_publication_and_conversion();
#ifdef EQUILIBRIUM_HAVE_NETCDF
        for (bool tracing : {false, true}) {
            test_netcdf_output(false, tracing);
            test_netcdf_output(true, tracing);
            if (tracing) {
                test_netcdf_output(false, tracing, true);
                test_netcdf_output(true, tracing, true);
                test_netcdf_output(false, tracing, true, true);
                test_netcdf_output(true, tracing, true, true);
                test_netcdf_output(false, tracing, true, true, true);
                test_netcdf_output(true, tracing, true, true, true);
            }
        }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
        for (bool tracing : {false, true}) {
            test_hdf5_output(false, tracing);
            test_hdf5_output(true, tracing);
            if (tracing) {
                test_hdf5_output(false, tracing, true);
                test_hdf5_output(true, tracing, true);
                test_hdf5_output(false, tracing, true, true);
                test_hdf5_output(true, tracing, true, true);
                test_hdf5_output(false, tracing, true, true, true);
                test_hdf5_output(true, tracing, true, true, true);
            }
        }
#endif
    } catch (const std::exception& error) {
        std::cerr << "boozer_binary_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "boozer_binary_test: PASS\n";
    return 0;
}
