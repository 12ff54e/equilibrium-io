#pragma once

#include "equilibrium/input_params.hpp"
#include "equilibrium/result.hpp"

#include <cstddef>
#include <filesystem>
#include <iosfwd>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace equilibrium {

struct LcfsPoint {
    double r = 0.0;
    double z = 0.0;
};

struct GfileLcfs {
    std::vector<LcfsPoint> points;  // Original polygon, metres.
    LcfsPoint magnetic_axis;
    int cocos = 1;  // Resolved source convention; coordinates stay in metres.
};

struct LcfsShape {
    double r_center = 0.0;  // (R_max + R_min)/2, not the magnetic-axis R.
    double minor_radius = 0.0;
    double z_reference = 0.0;
    double z_upper = 0.0;
    double z_lower = 0.0;
    double r_upper = 0.0;  // R at Z_max; midpoint of a flat top if tied.
    double r_lower = 0.0;
    double elongation_upper = 0.0;
    double elongation_lower = 0.0;
    double triangularity_upper = 0.0;
    double triangularity_lower = 0.0;
};

// Polygon extrema, without a spline fit. kappa = |Z_extreme-Z_reference|/a;
// delta = (R_center-R_at_extreme)/a. Positive delta points toward smaller R.
// The GfileLcfs overload uses magnetic_axis.z to distinguish the two halves.
LcfsShape calculate_lcfs_shape(std::span<const LcfsPoint> points,
                               double z_reference);
LcfsShape calculate_lcfs_shape(const GfileLcfs& lcfs);

// COCOS override: 1-8 or 11-18; otherwise use libmeq header/sign inference.
// Geometry reading does not require nonzero q/F or a star-shaped boundary.
// All functions except write_gfile_input throw on invalid input or I/O errors.
GfileLcfs read_gfile_lcfs(std::istream& input,
                          std::optional<int> cocos = std::nullopt);
GfileLcfs read_gfile_lcfs(const std::filesystem::path& path,
                          std::optional<int> cocos = std::nullopt);

struct GfileInputSettings {
    std::optional<int> cocos;
    int mpol = 12;   // cuMES supports 2-16.
    int ntheta = 0;  // Resolve to 2*mpol+6; otherwise even in [2*mpol+6,256].
    std::size_t boundary_points = 512;
    std::size_t profile_points = 129;
    double psi_ratio = 1.0;  // 1 preserves the supplied LCFS polygon.
    double boundary_tolerance = 1.0e-3;  // Maximum fit error / minor radius.
    double poloidal_flux_tolerance = 1.0e-3;  // Relative integrated error.
    // Native-double limits: ns in [3,512], max_iter>=1, ftol>=1e-16.
    InputStage stage{65, 10000, 1.0e-10};
};

struct GfileInput {
    InputParams input;  // Resolved fixed-boundary, fixed-iota cuMES input.
    GfileLcfs lcfs;     // Original LCFS, even when psi_ratio < 1.
    double psi_ratio = 1.0;
    double boundary_max_error = 0.0;  // Metres, relative to resampled contour.
    double poloidal_flux_relative_error = 0.0;
};

// Preserve physical field/current directions through COCOS 5. Pressure (Pa)
// and signed iota are lists on uniform normalized toroidal flux; aphi={1}.
// Conversion requires a boundary star-shaped about the magnetic axis and
// nonzero, sign-consistent q/F. A crop must bracket a closed flux contour.
GfileInput gfile_to_cumes_input(std::istream& input,
                                const GfileInputSettings& settings = {});
GfileInput gfile_to_cumes_input(const std::filesystem::path& path,
                                const GfileInputSettings& settings = {});

// Serialize the generated input as flat cuMES configuration JSON, not a
// scientific-result file. write_gfile_input publishes it atomically.
std::string gfile_input_to_json(const GfileInput& converted);
Status write_gfile_input(const std::string& path, const GfileInput& converted);

}  // namespace equilibrium
