#include "boozer_output_internal.hpp"
#include "equilibrium/writer_helpers.hpp"

#include <filesystem>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <netcdf.h>

namespace equilibrium::boozer_detail {
namespace {

void check_netcdf(int status, std::string_view operation) {
    if (status != NC_NOERR) {
        throw std::runtime_error("NetCDF " + std::string(operation) + ": " +
                                 nc_strerror(status));
    }
}

std::vector<double> periodic_angles(int count) {
    std::vector<double> values(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) {
        values[static_cast<std::size_t>(index)] = 2.0 * std::numbers::pi *
                                                  static_cast<double>(index) /
                                                  static_cast<double>(count);
    }
    return values;
}

}  // namespace

void write_boozer_netcdf(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    const RealBoozerSpectrum spectrum = make_real_spectrum(result);
    const auto temporary = reserve_output_path(path);
    std::error_code ignored;

    int file = -1;
    try {
        check_netcdf(nc_create(temporary.string().c_str(), NC_CLOBBER, &file),
                     "create");
        const int surfaces = result.source_ns - result.grid.first_surface;
        const int modes = static_cast<int>(spectrum.m.size());

        int surface_dim = -1;
        int theta_dim = -1;
        int alpha_dim = -1;
        int mode_dim = -1;
        check_netcdf(
            nc_def_dim(file, "surface", static_cast<std::size_t>(surfaces),
                       &surface_dim),
            "define surface dimension");
        check_netcdf(nc_def_dim(file, "theta_b",
                                static_cast<std::size_t>(result.grid.ntheta),
                                &theta_dim),
                     "define theta_b dimension");
        check_netcdf(
            nc_def_dim(file, "alpha",
                       static_cast<std::size_t>(result.grid.nzeta), &alpha_dim),
            "define alpha dimension");
        check_netcdf(nc_def_dim(file, "mode", static_cast<std::size_t>(modes),
                                &mode_dim),
                     "define mode dimension");

        const int surface_dims[1] = {surface_dim};
        const int theta_dims[1] = {theta_dim};
        const int alpha_dims[1] = {alpha_dim};
        const int mode_dims[1] = {mode_dim};
        const int spectrum_dims[2] = {surface_dim, mode_dim};
        const int grid_dims[3] = {surface_dim, alpha_dim, theta_dim};
        int s_var = -1;
        int iota_var = -1;
        int b2j00_var = -1;
        int theta_var = -1;
        int alpha_var = -1;
        int m_var = -1;
        int n_var = -1;
        int rmncc_var = -1;
        int rmnss_var = -1;
        int zmnsc_var = -1;
        int zmncs_var = -1;
        int numnsc_var = -1;
        int numncs_var = -1;
        int rmnsc_var = -1;
        int rmncs_var = -1;
        int zmncc_var = -1;
        int zmnss_var = -1;
        int numncc_var = -1;
        int numnss_var = -1;
        int b_var = -1;
        int sqrtg_var = -1;
        check_netcdf(nc_def_var(file, "s", NC_DOUBLE, 1, surface_dims, &s_var),
                     "define s");
        check_netcdf(
            nc_def_var(file, "iota", NC_DOUBLE, 1, surface_dims, &iota_var),
            "define iota");
        check_netcdf(
            nc_def_var(file, "b2j00", NC_DOUBLE, 1, surface_dims, &b2j00_var),
            "define b2j00");
        check_netcdf(
            nc_def_var(file, "theta_b", NC_DOUBLE, 1, theta_dims, &theta_var),
            "define theta_b");
        check_netcdf(
            nc_def_var(file, "alpha", NC_DOUBLE, 1, alpha_dims, &alpha_var),
            "define alpha");
        check_netcdf(nc_def_var(file, "mode_m", NC_INT, 1, mode_dims, &m_var),
                     "define mode_m");
        check_netcdf(nc_def_var(file, "mode_n", NC_INT, 1, mode_dims, &n_var),
                     "define mode_n");
        check_netcdf(
            nc_def_var(file, "rmncc", NC_DOUBLE, 2, spectrum_dims, &rmncc_var),
            "define rmncc");
        check_netcdf(
            nc_def_var(file, "rmnss", NC_DOUBLE, 2, spectrum_dims, &rmnss_var),
            "define rmnss");
        check_netcdf(
            nc_def_var(file, "zmnsc", NC_DOUBLE, 2, spectrum_dims, &zmnsc_var),
            "define zmnsc");
        check_netcdf(
            nc_def_var(file, "zmncs", NC_DOUBLE, 2, spectrum_dims, &zmncs_var),
            "define zmncs");
        check_netcdf(nc_def_var(file, "numnsc", NC_DOUBLE, 2, spectrum_dims,
                                &numnsc_var),
                     "define numnsc");
        check_netcdf(nc_def_var(file, "numncs", NC_DOUBLE, 2, spectrum_dims,
                                &numncs_var),
                     "define numncs");
        if (result.lasym) {
            check_netcdf(nc_def_var(file, "rmnsc", NC_DOUBLE, 2, spectrum_dims,
                                    &rmnsc_var),
                         "define rmnsc");
            check_netcdf(nc_def_var(file, "rmncs", NC_DOUBLE, 2, spectrum_dims,
                                    &rmncs_var),
                         "define rmncs");
            check_netcdf(nc_def_var(file, "zmncc", NC_DOUBLE, 2, spectrum_dims,
                                    &zmncc_var),
                         "define zmncc");
            check_netcdf(nc_def_var(file, "zmnss", NC_DOUBLE, 2, spectrum_dims,
                                    &zmnss_var),
                         "define zmnss");
            check_netcdf(nc_def_var(file, "numncc", NC_DOUBLE, 2, spectrum_dims,
                                    &numncc_var),
                         "define numncc");
            check_netcdf(nc_def_var(file, "numnss", NC_DOUBLE, 2, spectrum_dims,
                                    &numnss_var),
                         "define numnss");
        }
        check_netcdf(nc_def_var(file, "B", NC_DOUBLE, 3, grid_dims, &b_var),
                     "define B");
        check_netcdf(
            nc_def_var(file, "sqrt_g_b", NC_DOUBLE, 3, grid_dims, &sqrtg_var),
            "define sqrt_g_b");

        const auto put_text = [&](const char* name, std::string_view value) {
            check_netcdf(nc_put_att_text(file, NC_GLOBAL, name, value.size(),
                                         value.data()),
                         std::string("write attribute ") + name);
        };
        const auto put_int = [&](const char* name, int value) {
            check_netcdf(
                nc_put_att_int(file, NC_GLOBAL, name, NC_INT, 1, &value),
                std::string("write attribute ") + name);
        };
        const auto put_double = [&](const char* name, double value) {
            check_netcdf(
                nc_put_att_double(file, NC_GLOBAL, name, NC_DOUBLE, 1, &value),
                std::string("write attribute ") + name);
        };
        const bool has_flux_type =
            result.tracing && result.tracing->pflux_type == "list";
        const bool has_magnetic = result.tracing && result.tracing->magnetic;
        put_text("schema", has_flux_type     ? FLUX_BOOZER_SCHEMA
                           : result.profiles ? PROFILE_BOOZER_SCHEMA
                           : has_magnetic    ? MAGNETIC_BOOZER_SCHEMA
                           : result.tracing  ? TRACING_BOOZER_SCHEMA
                           : result.lasym    ? ASYMMETRIC_BOOZER_SCHEMA
                                             : BOOZER_SCHEMA);
        if (result.tracing || result.lasym)
            put_int("lasym", result.lasym ? 1 : 0);
        put_text("coordinate_convention", COORDINATE_CONVENTION);
        put_text("fourier_convention", result.lasym
                                           ? ASYMMETRIC_FOURIER_CONVENTION
                                           : FOURIER_CONVENTION);
        put_text("source_path", source_path.string());
        put_int("source_format_version", result.source_format_version);
        put_int("source_ns", result.source_ns);
        put_int("source_ntheta", result.source_ntheta);
        put_int("source_nzeta", result.source_nzeta);
        put_int("source_mpol", result.source_mpol);
        put_int("source_ntor", result.source_ntor);
        put_int("nfp", result.nfp);
        put_int("first_surface", result.grid.first_surface);
        put_int("ntheta", result.grid.ntheta);
        put_int("nzeta", result.grid.nzeta);
        put_int("mmax", result.spectrum.mmax);
        put_int("nmax", result.spectrum.nmax);
        put_int("radial_order", static_cast<int>(result.radial_order));
        put_double("resonance_tolerance", result.resonance_tolerance);

        std::vector<std::pair<int, std::span<const double>>> tracing_variables;
        if (result.tracing) {
            const auto& tracing = *result.tracing;
            put_text("tracing_convention", TRACING_CONVENTION);
            put_double("phiedge", tracing.phiedge);
            const auto dimension = [&](const char* name, std::size_t count) {
                int id = -1;
                check_netcdf(nc_def_dim(file, name, count, &id), name);
                return id;
            };
            const int full_dims[] = {dimension(
                "source_surface", static_cast<std::size_t>(result.source_ns))};
            const int aphi_dims[] = {
                dimension("aphi_coefficient", tracing.aphi.size())};
            const int axis_dims[] = {
                dimension("axis_mode", tracing.axis_r_cos.size())};
            const int half_dims[] = {
                dimension("axis_half_surface", tracing.axis_half_s.size())};
            const int native_dims[] = {
                half_dims[0],
                dimension("source_zeta",
                          static_cast<std::size_t>(result.source_nzeta)),
                dimension("source_theta",
                          static_cast<std::size_t>(result.source_ntheta))};
            const auto define = [&](const char* name, std::span<const int> dims,
                                    std::span<const double> values) {
                int id = -1;
                check_netcdf(
                    nc_def_var(file, name, NC_DOUBLE,
                               static_cast<int>(dims.size()), dims.data(), &id),
                    name);
                tracing_variables.emplace_back(id, values);
            };
            define("aphi", aphi_dims, tracing.aphi);
            put_text("pflux_type", tracing.pflux_type);
            define("toroidal_flux_derivative", full_dims,
                   tracing.toroidal_flux_derivative);
            define("poloidal_flux_derivative", full_dims,
                   tracing.poloidal_flux_derivative);
            define("axis_r_cos", axis_dims, tracing.axis_r_cos);
            define("axis_r_sin", axis_dims, tracing.axis_r_sin);
            define("axis_z_cos", axis_dims, tracing.axis_z_cos);
            define("axis_z_sin", axis_dims, tracing.axis_z_sin);
            define("axis_half_s", half_dims, tracing.axis_half_s);
            for (std::size_t field = 0; field < tracing.axis_half_fields.size();
                 ++field)
                define(AXIS_HALF_FIELD_NAMES[field], native_dims,
                       tracing.axis_half_fields[field]);
            if (tracing.magnetic) {
                const auto& magnetic = *tracing.magnetic;
                put_text("magnetic_convention", MAGNETIC_CONVENTION);
                for (std::size_t field = 0; field < magnetic.fields.size();
                     ++field)
                    define(MAGNETIC_FIELD_NAMES[field], grid_dims,
                           magnetic.fields[field]);
                define("source_theta_shift", grid_dims,
                       result.grid.source_theta_shift);
                define("nu_grid", grid_dims, result.grid.nu);
            }
        }

        if (result.profiles) {
            const auto& profiles = *result.profiles;
            put_text("profile_convention", PROFILE_CONVENTION);
            const auto define = [&](const char* name, const char* units,
                                    std::span<const double> values) {
                int id = -1;
                check_netcdf(
                    nc_def_var(file, name, NC_DOUBLE, 1, surface_dims, &id),
                    name);
                check_netcdf(
                    nc_put_att_text(file, id, "units",
                                    std::string_view(units).size(), units),
                    "write profile units");
                tracing_variables.emplace_back(id, values);
            };
            define("pressure", "Pa", profiles.pressure);
            define("current_i", "T m", profiles.current_i);
            define("current_g", "T m", profiles.current_g);
        }

        check_netcdf(nc_enddef(file), "end definitions");
        for (const auto& [id, values] : tracing_variables)
            check_netcdf(nc_put_var_double(file, id, values.data()),
                         "write tracing array");
        const auto theta_b = periodic_angles(result.grid.ntheta);
        const auto alpha = periodic_angles(result.grid.nzeta);
        check_netcdf(nc_put_var_double(file, s_var, result.s.data()),
                     "write s");
        check_netcdf(nc_put_var_double(file, iota_var, result.iota.data()),
                     "write iota");
        check_netcdf(
            nc_put_var_double(file, b2j00_var, result.grid.b2j00.data()),
            "write b2j00");
        check_netcdf(nc_put_var_double(file, theta_var, theta_b.data()),
                     "write theta_b");
        check_netcdf(nc_put_var_double(file, alpha_var, alpha.data()),
                     "write alpha");
        check_netcdf(nc_put_var_int(file, m_var, spectrum.m.data()),
                     "write mode_m");
        check_netcdf(nc_put_var_int(file, n_var, spectrum.n.data()),
                     "write mode_n");
        check_netcdf(nc_put_var_double(file, rmncc_var, spectrum.rmncc.data()),
                     "write rmncc");
        check_netcdf(nc_put_var_double(file, rmnss_var, spectrum.rmnss.data()),
                     "write rmnss");
        check_netcdf(nc_put_var_double(file, zmnsc_var, spectrum.zmnsc.data()),
                     "write zmnsc");
        check_netcdf(nc_put_var_double(file, zmncs_var, spectrum.zmncs.data()),
                     "write zmncs");
        check_netcdf(
            nc_put_var_double(file, numnsc_var, spectrum.numnsc.data()),
            "write numnsc");
        check_netcdf(
            nc_put_var_double(file, numncs_var, spectrum.numncs.data()),
            "write numncs");
        if (result.lasym) {
            check_netcdf(
                nc_put_var_double(file, rmnsc_var, spectrum.rmnsc.data()),
                "write rmnsc");
            check_netcdf(
                nc_put_var_double(file, rmncs_var, spectrum.rmncs.data()),
                "write rmncs");
            check_netcdf(
                nc_put_var_double(file, zmncc_var, spectrum.zmncc.data()),
                "write zmncc");
            check_netcdf(
                nc_put_var_double(file, zmnss_var, spectrum.zmnss.data()),
                "write zmnss");
            check_netcdf(
                nc_put_var_double(file, numncc_var, spectrum.numncc.data()),
                "write numncc");
            check_netcdf(
                nc_put_var_double(file, numnss_var, spectrum.numnss.data()),
                "write numnss");
        }
        check_netcdf(nc_put_var_double(file, b_var, result.grid.b.data()),
                     "write B");
        check_netcdf(
            nc_put_var_double(file, sqrtg_var, result.grid.sqrtg_b.data()),
            "write sqrt_g_b");
        check_netcdf(nc_close(file), "close");
        file = -1;
        const auto error =
            io_detail::publish_library_file(temporary.string(), path.string());
        if (!error.empty()) throw std::runtime_error(error);
    } catch (...) {
        if (file >= 0) nc_close(file);
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

}  // namespace equilibrium::boozer_detail
