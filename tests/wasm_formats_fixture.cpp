#include "equilibrium/boozer_binary.hpp"
#include "equilibrium/writer.hpp"
#include "input_fixture.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>

#include <emscripten.h>
#include <netcdf.h>

namespace {

void equilibrium_fixture(int kind) {
    auto params = equilibrium::test::input_fixture();
    params.mpol = 3;
    params.ntor = 1;
    params.nfp = 2;
    params.ntheta = 12;
    params.nzeta = 6;
    params.lasym = kind == 1 || kind == 3 || kind == 10;
    params.stages = {{3, 10, 1e-5}, {5, 20, 1e-9}};
    params.rbc_m = {0, 1, 2, 2};
    params.rbc_n = {0, 0, -1, -1};
    params.rbc_value = {4.0, 1.0, 0.01, -0.003};
    params.zbs_m = {1, 1};
    params.zbs_n = {0, 1};
    params.zbs_value = {1.0, 0.0};
    params.rbcc = {4, 0, 1, 0, 0, 0.01 - 0.003};
    params.rbss = {0, 0, 0, 0, 0, -0.01 + 0.003};
    params.zbsc = {0, 0, 1, 0, 0, 0};
    params.zbcs.assign(6, 0.0);
    params.pmass_type = "two_power";
    params.am = {0.2, 1.0, 2.0};
    params.ai = {0.4, 0.1};
    params.aphi = {1.0, 0.2};
    params.raxis_c = {4.0, 0.01};
    params.zaxis_s = {0.0, 0.0};
    if (kind == 9 || kind == 10) {
        params.pmass_type = params.piota_type = params.pcurr_type =
            params.pflux_type = "list";
        params.am = {8.0, 2.0, 6.0};
        params.ai = {-0.2, 0.4};
        params.ac = {-1, 2};
        params.aphi = {1, 2};
    }
    if (params.lasym) {
        params.raxis_s = {0.0, 0.002};
        params.zaxis_c = {0.03, -0.001};
        params.rbsc.assign(6, 0.0);
        params.rbcs.assign(6, 0.0);
        params.zbcc.assign(6, 0.0);
        params.zbss.assign(6, 0.0);
    }
    if (kind == 1) {
        params.rbs_m = {1, 1, 2};
        params.rbs_n = {-1, -1, 1};
        params.rbs_value = {0.02, -0.003, 0.0};
        params.zbc_m = {0, 2};
        params.zbc_n = {0, 1};
        params.zbc_value = {0.03, 0.01};
        params.rbsc[3] = params.rbcs[3] = 0.02 - 0.003;
        params.zbcc[0] = 0.03;
        params.zbcc[5] = params.zbss[5] = 0.01;
    }
    if (kind == 2) {
        params.lfreeb = true;
        params.nvacskip = 3;
        params.coils_file = "/inputs/coils.test";
        params.extcur = {1.25, -2.0};
        params.embedded_makegrid_parameters =
            equilibrium::MakegridParametersSpec{false, false, 2,   1.0, 7.0,
                                                5,     -3.0,  3.0, 5,   4};
    }
    equilibrium::EquilibriumSnapshot state;
    state.ns = 5;
    state.mnmax = 6;
    state.families.resize(params.lasym ? 12 : 6);
    for (int c = 0; c < state.components(); ++c) {
        auto& family = state.families[static_cast<std::size_t>(c)];
        family.resize(state.family_size());
        for (std::size_t i = 0; i < family.size(); ++i)
            family[i] = 0.25 * c + 0.001 * static_cast<double>(i);
    }
    state.ntheta = params.ntheta;
    state.nzeta = params.nzeta;
    for (std::size_t c = 0; c < state.half_fields.size(); ++c) {
        auto& field = state.half_fields[c];
        field.resize(state.half_field_size());
        for (std::size_t i = 0; i < field.size(); ++i)
            field[i] = -2.0 - 0.25 * c - 0.0001 * static_cast<double>(i);
    }
    for (std::size_t c = 0; c < state.full_fields.size(); ++c)
        state.full_fields[c].assign(state.full_field_size(), 0.125 * (c + 1));
    equilibrium::RunReport report;
    report.status = equilibrium::RunStatus::CONVERGED;
    report.total_effective_iterations = 17;
    report.stages = {{3, 6, true, {1e-6, 2e-6, 3e-6}, {{1}, {4}}},
                     {5, 11, true, {1e-11, 2e-11, 3e-11}, {{3}}}};
    report.build = {"original-revision", true,
                    "Release",           "float",
                    "mixed-float",       "original computation flags"};
    report.input = {"/inputs/field-µ.json", "original-source-hash"};
    report.runtime = {"Test GPU", "test-driver", "WebGPU", "original-toolkit"};
    report.input_params = params;
    if (kind == 6) report.input_params.rbcc[0] += 1.0;
    const auto status = equilibrium::make_binary_writer()->write_atomic(
        equilibrium::as_view(state), report,
        {equilibrium::OutputFormat::BINARY, "/fixture.bin"});
    if (!status) throw std::runtime_error(status.error());
}

void boozer_fixture(bool lasym, bool with_profiles, bool flux_list) {
    using namespace equilibrium;
    BoozerResult result;
    result.lasym = lasym;
    result.source_format_version = lasym ? 10 : 8;
    result.source_ns = 3;
    result.source_ntheta = 8;
    result.source_nzeta = 6;
    result.source_mpol = 3;
    result.source_ntor = 1;
    result.nfp = 2;
    result.radial_order = RadialInterpolationOrder::TWO_POINT;
    result.s = {0.5, 1.0};
    result.iota = {0.4, 0.5};
    result.grid.source_ns = result.spectrum.source_ns = 3;
    result.grid.first_surface = result.spectrum.first_surface = 1;
    result.grid.ntheta = result.grid.nzeta = 4;
    result.grid.b.resize(32);
    result.grid.sqrtg_b.resize(32);
    for (std::size_t i = 0; i < 32; ++i) {
        result.grid.b[i] = static_cast<double>(i + 1);
        result.grid.sqrtg_b[i] = -static_cast<double>(i + 1);
    }
    result.grid.b2j00 = {-10, -11};
    result.grid.source_theta_shift.assign(32, 4.125);
    result.grid.nu.assign(32, 0.071);
    result.spectrum.mmax = result.spectrum.nmax = 1;
    for (int n = -1; n <= 1; ++n)
        for (int m = -1; m <= 1; ++m) {
            result.spectrum.m.push_back(m);
            result.spectrum.n.push_back(n);
        }
    for (int surface = 1; surface <= 2; ++surface)
        for (int n = -1; n <= 1; ++n)
            for (int m = -1; m <= 1; ++m) {
                // Near-cancelling parity pairs also exercise re-encoding's
                // absolute error bound when one output coefficient is tiny.
                result.spectrum.r.emplace_back(
                    surface * (1.23456789012345 +
                               0.123456789012345 * std::abs(m) + 1e-16 * m * n),
                    lasym ? surface * 0.02 * (m + 2 * n) : 0.0);
                result.spectrum.z.emplace_back(
                    lasym ? surface * 0.125 * (1 + std::abs(n)) : 0.0,
                    surface * 0.3 * (m + 2 * n));
                result.spectrum.nu.emplace_back(
                    lasym ? surface * 0.01 * (1 + std::abs(m)) : 0.0,
                    surface * 0.03 * (m - 3 * n));
            }
    auto& tracing = result.tracing.emplace();
    tracing.phiedge = 7.5;
    tracing.aphi = {1, 0.2};
    tracing.toroidal_flux_derivative = {-2, -3, -4};
    tracing.poloidal_flux_derivative = {-0.6, -1.2, -2};
    tracing.axis_r_cos = {3, 0.1};
    tracing.axis_r_sin = {0, lasym ? 0.02 : 0};
    tracing.axis_z_cos = {lasym ? 0.12 : 0, 0};
    tracing.axis_z_sin = {0, 0.03};
    tracing.axis_half_s = {0.25, 0.75};
    for (auto& field : tracing.axis_half_fields) field.assign(2 * 8 * 6, -2.5);
    auto& magnetic = tracing.magnetic.emplace();
    for (auto& field : magnetic.fields) field.assign(32, -1.5);
    if (with_profiles) {
        auto& profiles = result.profiles.emplace();
        profiles.pressure = {125000.0, 0.0};
        profiles.current_i = {1.25, 2.5};
        profiles.current_g = {-5.0, -7.0};
    }
    if (flux_list) {
        result.source_format_version = 11;
        tracing.pflux_type = "list";
        tracing.aphi = {1, 2, 1};
    }
    write_boozer_binary("/fixture.bin", result, "original-equilibrium.bin");
}

}  // namespace

extern "C" EMSCRIPTEN_KEEPALIVE int equilibrium_test_format_fixture(int kind) {
    try {
        if (kind < 4 || kind == 6 || kind == 9 || kind == 10)
            equilibrium_fixture(kind);
        else
            boozer_fixture(kind == 5 || kind == 8 || kind == 12, kind >= 7,
                           kind >= 11);
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}

extern "C" EMSCRIPTEN_KEEPALIVE int equilibrium_test_netcdf4() {
    int file = -1, dim = -1, variable = -1;
    const double values[3] = {1.25, -2.5, 3.75};
    double decoded[3]{};
    if (nc_create("/compressed.nc", NC_NETCDF4 | NC_CLOBBER, &file) ||
        nc_def_dim(file, "sample", 3, &dim) ||
        nc_def_var(file, "sample", NC_DOUBLE, 1, &dim, &variable) ||
        nc_def_var_deflate(file, variable, 1, 1, 4) || nc_enddef(file) ||
        nc_put_var_double(file, variable, values) || nc_close(file) ||
        nc_open("/compressed.nc", NC_NOWRITE, &file) ||
        nc_inq_varid(file, "sample", &variable) ||
        nc_get_var_double(file, variable, decoded) || nc_close(file))
        return 1;
    for (int i = 0; i < 3; ++i)
        if (decoded[i] != values[i]) return 1;
    std::remove("/compressed.nc");
    return 0;
}

// Produce a compressed NetCDF-4 container through the dependency's C API,
// independently of the cuMES writer. Retain exact attribute types and bytes.
extern "C" EMSCRIPTEN_KEEPALIVE int equilibrium_test_netcdf4_copy() {
    int source = -1, output = -1;
    const auto check = [](int code) {
        if (code) throw std::runtime_error(nc_strerror(code));
    };
    try {
        check(nc_open("/fixture.nc", NC_NOWRITE, &source));
        check(nc_create("/fixture4.nc", NC_CLOBBER | NC_NETCDF4, &output));
        int dimensions = 0, variables = 0, attributes = 0, unlimited = -1;
        check(nc_inq(source, &dimensions, &variables, &attributes, &unlimited));
        std::vector<int> dims(static_cast<std::size_t>(dimensions));
        for (int d = 0; d < dimensions; ++d) {
            char name[NC_MAX_NAME + 1]{};
            std::size_t count = 0;
            check(nc_inq_dim(source, d, name, &count));
            check(nc_def_dim(output, name,
                             d == unlimited ? NC_UNLIMITED : count,
                             &dims[static_cast<std::size_t>(d)]));
        }
        for (int a = 0; a < attributes; ++a) {
            char name[NC_MAX_NAME + 1]{};
            check(nc_inq_attname(source, NC_GLOBAL, a, name));
            check(nc_copy_att(source, NC_GLOBAL, name, output, NC_GLOBAL));
        }
        for (int v = 0; v < variables; ++v) {
            char name[NC_MAX_NAME + 1]{};
            int rank = 0, id = -1;
            nc_type type;
            check(nc_inq_var(source, v, name, &type, &rank, nullptr, nullptr));
            std::vector<int> shape(static_cast<std::size_t>(rank));
            check(nc_inq_vardimid(source, v, shape.data()));
            std::size_t count = 1, width = 0;
            check(nc_inq_type(source, type, nullptr, &width));
            for (auto& dim : shape) {
                std::size_t length = 0;
                check(nc_inq_dimlen(source, dim, &length));
                count *= length;
                dim = dims[static_cast<std::size_t>(dim)];
            }
            check(nc_def_var(output, name, type, rank, shape.data(), &id));
            if (rank) check(nc_def_var_deflate(output, id, 1, 1, 4));
            std::vector<unsigned char> values(count * width);
            check(nc_get_var(source, v, values.data()));
            check(nc_put_var(output, id, values.data()));
        }
        check(nc_close(source));
        source = -1;
        check(nc_close(output));
        output = -1;
        return 0;
    } catch (const std::exception& error) {
        if (source >= 0) nc_close(source);
        if (output >= 0) nc_close(output);
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
