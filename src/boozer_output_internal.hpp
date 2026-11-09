#ifndef EQUILIBRIUM_SRC_BOOZER_OUTPUT_INTERNAL_HPP_
#define EQUILIBRIUM_SRC_BOOZER_OUTPUT_INTERNAL_HPP_

#include "equilibrium/boozer.hpp"
#include "equilibrium/writer_helpers.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace equilibrium::boozer_detail {

// Reserve the path before a library opens it: even a failed HDF5/NetCDF
// create can leave a partial file, which the caller then owns and removes.
inline std::filesystem::path reserve_output_path(
    const std::filesystem::path& path) {
    const auto temporary = io_detail::temp_path_for(path.string());
    auto* file = std::fopen(temporary.c_str(), "wbx");
    if (!file) throw std::runtime_error("cannot reserve Boozer output");
    if (std::fclose(file) != 0) {
        std::remove(temporary.c_str());
        throw std::runtime_error("cannot close reserved Boozer output");
    }
    return temporary;
}

inline constexpr std::string_view BOOZER_SCHEMA =
    "magnetic-coordinate-boozer-v8";
inline constexpr std::string_view PROFILE_CONVENTION =
    "surface-profiles-v1: pressure in Pa; B=beta*grad(s)+I*grad(theta_b)"
    "+G*grad(zeta_b); current_i=I and current_g=G in T m with physical "
    "zeta_b; all profiles share the non-axis s/iota grid";
inline constexpr std::string_view MAGNETIC_CONVENTION =
    "boozer-magnetic-v1: boozer fields use (s,theta_b,zeta_b) with physical "
    "zeta_b; sampled on the mixed grid with theta=theta_b+source_theta_shift "
    "and zeta_b=alpha/nfp+nu_grid; shifts are unwrapped and untruncated";
inline constexpr std::array<const char*, EquilibriumSnapshot::HALF_FIELD_COUNT>
    MAGNETIC_FIELD_NAMES = {"boozer_sqrtg", "boozer_bsups", "boozer_bsupu",
                            "boozer_bsupv", "boozer_bsubs", "boozer_bsubu",
                            "boozer_bsubv"};
inline constexpr std::string_view TRACING_CONVENTION =
    "tracing-v1: signed flux/(2*pi) derivatives in source s; axis m=0 "
    "geometry in alpha=nfp*zeta; axis_half fields are unmodified native "
    "half-grid samples in (s,theta,zeta), not Boozer axis values";
inline constexpr std::array<const char*, EquilibriumSnapshot::HALF_FIELD_COUNT>
    AXIS_HALF_FIELD_NAMES = {"axis_half_sqrtg", "axis_half_bsups",
                             "axis_half_bsupu", "axis_half_bsupv",
                             "axis_half_bsubs", "axis_half_bsubu",
                             "axis_half_bsubv"};
inline constexpr std::string_view COORDINATE_CONVENTION =
    "mixed-grid-v2: theta_b uniform; alpha=nfp*zeta is the unchanged source "
    "field-period angle; alpha_b=alpha+nfp*nu; nu is in physical toroidal "
    "radians";
inline constexpr std::string_view FOURIER_CONVENTION =
    "real-parity-v3: f=sum[cc*cos(m*theta_b)*cos(n*alpha) + "
    "ss*sin(m*theta_b)*sin(n*alpha)] for even fields and "
    "f=sum[sc*sin(m*theta_b)*cos(n*alpha) + "
    "cs*cos(m*theta_b)*sin(n*alpha)] for odd fields; m,n are nonnegative "
    "and alpha is a field-period angle";
inline constexpr std::string_view ASYMMETRIC_FOURIER_CONVENTION =
    "real-parity-v4: f=sum[cc*cos(m*theta_b)*cos(n*alpha) + "
    "ss*sin(m*theta_b)*sin(n*alpha) + sc*sin(m*theta_b)*cos(n*alpha) + "
    "cs*cos(m*theta_b)*sin(n*alpha)]; m,n are nonnegative "
    "and alpha is a field-period angle; no stellarator symmetry is imposed";

using ::equilibrium::RealBoozerSpectrum;

void validate_result(const BoozerResult& result);
void validate_group(const BoozerResult& result, DataGroup group);
void validate_real_spectrum(const RealBoozerSpectrum& spectrum,
                            int surfaces,
                            bool geometry = true,
                            bool shift = true);
RealBoozerSpectrum make_real_spectrum(const BoozerResult& result);
MixedGridSpectrum restore_complex_spectrum(const RealBoozerSpectrum& spectrum,
                                           int source_ns,
                                           int first_surface,
                                           bool geometry = true,
                                           bool shift = true);

#ifdef EQUILIBRIUM_HAVE_NETCDF
void write_boozer_netcdf(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);
#endif

#ifdef EQUILIBRIUM_HAVE_HDF5
void write_boozer_hdf5(const std::filesystem::path& path,
                       const BoozerResult& result,
                       const std::filesystem::path& source_path);
#endif

}  // namespace equilibrium::boozer_detail

#endif  // EQUILIBRIUM_SRC_BOOZER_OUTPUT_INTERNAL_HPP_
