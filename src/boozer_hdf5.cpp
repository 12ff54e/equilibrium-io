#include "boozer_output_internal.hpp"
#include "hdf5_output.hpp"

#include <filesystem>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <hdf5.h>

namespace equilibrium::boozer_detail {
namespace {

void check_hdf5(herr_t status, std::string_view operation) {
    if (status < 0) {
        throw std::runtime_error("HDF5 " + std::string(operation) + " failed");
    }
}

using io_detail::Hdf5Handle;

void write_string_attribute(hid_t file,
                            const char* name,
                            const std::string& value) {
    Hdf5Handle<H5Tclose> type(H5Tcopy(H5T_C_S1));
    if (type.get() < 0)
        throw std::runtime_error("HDF5 copy string type failed");
    check_hdf5(H5Tset_size(type.get(), value.size() + 1), "set string width");
    Hdf5Handle<H5Sclose> space(H5Screate(H5S_SCALAR));
    if (space.get() < 0)
        throw std::runtime_error("HDF5 create scalar space failed");
    Hdf5Handle<H5Aclose> attribute(H5Acreate2(
        file, name, type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT));
    if (attribute.get() < 0)
        throw std::runtime_error("HDF5 create string attribute failed");
    check_hdf5(H5Awrite(attribute.get(), type.get(), value.c_str()),
               "write string attribute");
    check_hdf5(attribute.close(), "close string attribute");
}

template <typename T>
void write_scalar_attribute(hid_t file,
                            const char* name,
                            hid_t type,
                            const T& value) {
    Hdf5Handle<H5Sclose> space(H5Screate(H5S_SCALAR));
    if (space.get() < 0)
        throw std::runtime_error("HDF5 create scalar space failed");
    Hdf5Handle<H5Aclose> attribute(
        H5Acreate2(file, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT));
    if (attribute.get() < 0)
        throw std::runtime_error("HDF5 create scalar attribute failed");
    check_hdf5(H5Awrite(attribute.get(), type, &value),
               "write scalar attribute");
    check_hdf5(attribute.close(), "close scalar attribute");
}

void write_dataset(hid_t file,
                   const char* name,
                   hid_t type,
                   std::span<const hsize_t> dimensions,
                   const void* data,
                   const char* units = nullptr) {
    Hdf5Handle<H5Sclose> space(H5Screate_simple(
        static_cast<int>(dimensions.size()), dimensions.data(), nullptr));
    if (space.get() < 0)
        throw std::runtime_error("HDF5 create dataset space failed");
    Hdf5Handle<H5Dclose> dataset(H5Dcreate2(
        file, name, type, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    if (dataset.get() < 0)
        throw std::runtime_error("HDF5 create dataset failed");
    check_hdf5(
        H5Dwrite(dataset.get(), type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data),
        "write dataset");
    if (units) write_string_attribute(dataset.get(), "units", units);
    check_hdf5(dataset.close(), "close dataset");
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

void write_boozer_hdf5(const std::filesystem::path& path,
                       const BoozerResult& result,
                       const std::filesystem::path& source_path) {
    const RealBoozerSpectrum spectrum = make_real_spectrum(result);
    io_detail::Hdf5Output output(path.string());
    const hid_t file = output.get();
    if (file < 0) throw std::runtime_error("HDF5 create image failed");

    const bool has_flux_type =
        result.tracing && result.tracing->pflux_type == "list";
    const bool has_magnetic = result.tracing && result.tracing->magnetic;
    write_string_attribute(file, "schema",
                           std::string(has_flux_type     ? FLUX_BOOZER_SCHEMA
                                       : result.profiles ? PROFILE_BOOZER_SCHEMA
                                       : has_magnetic   ? MAGNETIC_BOOZER_SCHEMA
                                       : result.tracing ? TRACING_BOOZER_SCHEMA
                                       : result.lasym ? ASYMMETRIC_BOOZER_SCHEMA
                                                      : BOOZER_SCHEMA));
    write_string_attribute(file, "coordinate_convention",
                           std::string(COORDINATE_CONVENTION));
    write_string_attribute(
        file, "fourier_convention",
        std::string(result.lasym ? ASYMMETRIC_FOURIER_CONVENTION
                                 : FOURIER_CONVENTION));
    if (result.tracing || result.lasym)
        write_scalar_attribute(file, "lasym", H5T_NATIVE_INT,
                               result.lasym ? 1 : 0);
    write_string_attribute(file, "source_path", source_path.string());
    write_scalar_attribute(file, "source_format_version", H5T_NATIVE_INT,
                           result.source_format_version);
    write_scalar_attribute(file, "source_ns", H5T_NATIVE_INT, result.source_ns);
    write_scalar_attribute(file, "source_ntheta", H5T_NATIVE_INT,
                           result.source_ntheta);
    write_scalar_attribute(file, "source_nzeta", H5T_NATIVE_INT,
                           result.source_nzeta);
    write_scalar_attribute(file, "source_mpol", H5T_NATIVE_INT,
                           result.source_mpol);
    write_scalar_attribute(file, "source_ntor", H5T_NATIVE_INT,
                           result.source_ntor);
    write_scalar_attribute(file, "nfp", H5T_NATIVE_INT, result.nfp);
    write_scalar_attribute(file, "first_surface", H5T_NATIVE_INT,
                           result.grid.first_surface);
    write_scalar_attribute(file, "ntheta", H5T_NATIVE_INT, result.grid.ntheta);
    write_scalar_attribute(file, "nzeta", H5T_NATIVE_INT, result.grid.nzeta);
    write_scalar_attribute(file, "mmax", H5T_NATIVE_INT, result.spectrum.mmax);
    write_scalar_attribute(file, "nmax", H5T_NATIVE_INT, result.spectrum.nmax);
    const int radial_order = static_cast<int>(result.radial_order);
    write_scalar_attribute(file, "radial_order", H5T_NATIVE_INT, radial_order);
    write_scalar_attribute(file, "resonance_tolerance", H5T_NATIVE_DOUBLE,
                           result.resonance_tolerance);

    if (result.tracing) {
        const auto& tracing = *result.tracing;
        write_string_attribute(file, "tracing_convention",
                               std::string(TRACING_CONVENTION));
        write_scalar_attribute(file, "phiedge", H5T_NATIVE_DOUBLE,
                               tracing.phiedge);
        const hsize_t full_dims[] = {static_cast<hsize_t>(result.source_ns)};
        const hsize_t aphi_dims[] = {static_cast<hsize_t>(tracing.aphi.size())};
        const hsize_t axis_dims[] = {
            static_cast<hsize_t>(tracing.axis_r_cos.size())};
        const hsize_t half_dims[] = {
            static_cast<hsize_t>(tracing.axis_half_s.size())};
        const hsize_t native_dims[] = {
            half_dims[0], static_cast<hsize_t>(result.source_nzeta),
            static_cast<hsize_t>(result.source_ntheta)};
        write_dataset(file, "aphi", H5T_NATIVE_DOUBLE, aphi_dims,
                      tracing.aphi.data());
        write_string_attribute(file, "pflux_type", tracing.pflux_type);
        write_dataset(file, "toroidal_flux_derivative", H5T_NATIVE_DOUBLE,
                      full_dims, tracing.toroidal_flux_derivative.data());
        write_dataset(file, "poloidal_flux_derivative", H5T_NATIVE_DOUBLE,
                      full_dims, tracing.poloidal_flux_derivative.data());
        write_dataset(file, "axis_r_cos", H5T_NATIVE_DOUBLE, axis_dims,
                      tracing.axis_r_cos.data());
        write_dataset(file, "axis_r_sin", H5T_NATIVE_DOUBLE, axis_dims,
                      tracing.axis_r_sin.data());
        write_dataset(file, "axis_z_cos", H5T_NATIVE_DOUBLE, axis_dims,
                      tracing.axis_z_cos.data());
        write_dataset(file, "axis_z_sin", H5T_NATIVE_DOUBLE, axis_dims,
                      tracing.axis_z_sin.data());
        write_dataset(file, "axis_half_s", H5T_NATIVE_DOUBLE, half_dims,
                      tracing.axis_half_s.data());
        for (std::size_t field = 0; field < tracing.axis_half_fields.size();
             ++field)
            write_dataset(file, AXIS_HALF_FIELD_NAMES[field], H5T_NATIVE_DOUBLE,
                          native_dims, tracing.axis_half_fields[field].data());
        if (tracing.magnetic) {
            const auto& magnetic = *tracing.magnetic;
            write_string_attribute(file, "magnetic_convention",
                                   std::string(MAGNETIC_CONVENTION));
            const hsize_t map_dims[] = {
                static_cast<hsize_t>(result.s.size()),
                static_cast<hsize_t>(result.grid.nzeta),
                static_cast<hsize_t>(result.grid.ntheta)};
            for (std::size_t field = 0; field < magnetic.fields.size(); ++field)
                write_dataset(file, MAGNETIC_FIELD_NAMES[field],
                              H5T_NATIVE_DOUBLE, map_dims,
                              magnetic.fields[field].data());
            write_dataset(file, "source_theta_shift", H5T_NATIVE_DOUBLE,
                          map_dims, result.grid.source_theta_shift.data());
            write_dataset(file, "nu_grid", H5T_NATIVE_DOUBLE, map_dims,
                          result.grid.nu.data());
        }
    }

    const hsize_t surface_dims[1] = {static_cast<hsize_t>(result.s.size())};
    const hsize_t theta_dims[1] = {static_cast<hsize_t>(result.grid.ntheta)};
    const hsize_t alpha_dims[1] = {static_cast<hsize_t>(result.grid.nzeta)};
    const hsize_t mode_dims[1] = {static_cast<hsize_t>(spectrum.m.size())};
    const hsize_t spectrum_dims[2] = {static_cast<hsize_t>(result.s.size()),
                                      static_cast<hsize_t>(spectrum.m.size())};
    const hsize_t grid_dims[3] = {static_cast<hsize_t>(result.s.size()),
                                  static_cast<hsize_t>(result.grid.nzeta),
                                  static_cast<hsize_t>(result.grid.ntheta)};
    const auto theta_b = periodic_angles(result.grid.ntheta);
    const auto alpha = periodic_angles(result.grid.nzeta);
    write_dataset(file, "s", H5T_NATIVE_DOUBLE, surface_dims, result.s.data());
    write_dataset(file, "iota", H5T_NATIVE_DOUBLE, surface_dims,
                  result.iota.data());
    if (result.profiles) {
        const auto& profiles = *result.profiles;
        write_string_attribute(file, "profile_convention",
                               std::string(PROFILE_CONVENTION));
        write_dataset(file, "pressure", H5T_NATIVE_DOUBLE, surface_dims,
                      profiles.pressure.data(), "Pa");
        write_dataset(file, "current_i", H5T_NATIVE_DOUBLE, surface_dims,
                      profiles.current_i.data(), "T m");
        write_dataset(file, "current_g", H5T_NATIVE_DOUBLE, surface_dims,
                      profiles.current_g.data(), "T m");
    }
    write_dataset(file, "b2j00", H5T_NATIVE_DOUBLE, surface_dims,
                  result.grid.b2j00.data());
    write_dataset(file, "theta_b", H5T_NATIVE_DOUBLE, theta_dims,
                  theta_b.data());
    write_dataset(file, "alpha", H5T_NATIVE_DOUBLE, alpha_dims, alpha.data());
    write_dataset(file, "mode_m", H5T_NATIVE_INT, mode_dims, spectrum.m.data());
    write_dataset(file, "mode_n", H5T_NATIVE_INT, mode_dims, spectrum.n.data());
    write_dataset(file, "rmncc", H5T_NATIVE_DOUBLE, spectrum_dims,
                  spectrum.rmncc.data());
    write_dataset(file, "rmnss", H5T_NATIVE_DOUBLE, spectrum_dims,
                  spectrum.rmnss.data());
    write_dataset(file, "zmnsc", H5T_NATIVE_DOUBLE, spectrum_dims,
                  spectrum.zmnsc.data());
    write_dataset(file, "zmncs", H5T_NATIVE_DOUBLE, spectrum_dims,
                  spectrum.zmncs.data());
    write_dataset(file, "numnsc", H5T_NATIVE_DOUBLE, spectrum_dims,
                  spectrum.numnsc.data());
    write_dataset(file, "numncs", H5T_NATIVE_DOUBLE, spectrum_dims,
                  spectrum.numncs.data());
    if (result.lasym) {
        write_dataset(file, "rmnsc", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.rmnsc.data());
        write_dataset(file, "rmncs", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.rmncs.data());
        write_dataset(file, "zmncc", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.zmncc.data());
        write_dataset(file, "zmnss", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.zmnss.data());
        write_dataset(file, "numncc", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.numncc.data());
        write_dataset(file, "numnss", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.numnss.data());
    }
    write_dataset(file, "B", H5T_NATIVE_DOUBLE, grid_dims,
                  result.grid.b.data());
    write_dataset(file, "sqrt_g_b", H5T_NATIVE_DOUBLE, grid_dims,
                  result.grid.sqrtg_b.data());
    const auto error = output.publish();
    if (!error.empty()) throw std::runtime_error("HDF5: " + error);
}

}  // namespace equilibrium::boozer_detail
