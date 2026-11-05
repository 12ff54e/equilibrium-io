#pragma once

#include "equilibrium/equilibrium_snapshot.hpp"

#include <array>
#include <complex>
#include <optional>
#include <string>
#include <vector>

namespace equilibrium {

// Boozer records use the existing mixed grid: theta_b is uniform and
// alpha=nfp*zeta remains the source field-period angle. These are data records;
// coordinate transforms and interpolation operations belong to consumers.
enum class RadialInterpolationOrder {
    TWO_POINT = 2,
    FOUR_POINT = 4,
};

struct BoozerMixedGrid {
    int source_ns = 0;
    int first_surface = 1;
    int ntheta = 0;
    int nzeta = 0;

    std::vector<double> source_theta_p;
    std::vector<double> source_theta;
    // Smooth periodic displacement theta_native - theta_b, without wrapping.
    // Together with nu, this specifies the full source-to-Boozer map.
    std::vector<double> source_theta_shift;
    std::vector<double> r;
    std::vector<double> z;
    std::vector<double> nu;
    std::vector<double> b;
    std::vector<double> sqrtg_b;
    std::vector<double> b2j00;
};

struct MixedGridSpectrum {
    int source_ns = 0;
    int first_surface = 1;
    int mmax = 0;
    int nmax = 0;

    // Mode ordering is n-major, then m from -mmax through mmax.
    std::vector<int> m;
    std::vector<int> n;
    std::vector<std::complex<double>> r;
    std::vector<std::complex<double>> z;
    std::vector<std::complex<double>> nu;
};

struct BoozerMagneticData {
    // Signed Jacobian and six components in (s, theta_b, physical zeta_b).
    // Ordering matches EquilibriumSnapshot::HALF_FIELD_COUNT; sampling matches
    // BoozerMixedGrid: [non-axis surface][source alpha][theta_b].
    std::array<std::vector<double>, EquilibriumSnapshot::HALF_FIELD_COUNT>
        fields;
};

struct BoozerTracingData {
    // Original cuMES flux-profile inputs (SI), without a polarity change.
    double phiedge = 0.0;
    std::vector<double> aphi;
    std::string pflux_type = "power_series";

    // Signed <sqrt(g) B^zeta> and <sqrt(g) B^theta>, interpolated to all
    // source full-grid surfaces s=j/(source_ns-1), including the axis.
    // zeta is the physical toroidal angle; these are flux/(2*pi) derivatives.
    std::vector<double> toroidal_flux_derivative;
    std::vector<double> poloidal_flux_derivative;

    // Converged m=0 axis geometry, all n=0..source_ntor, independent of the
    // output truncation: R/Z=sum[cos[n]*cos(n*alpha)+sin[n]*sin(n*alpha)].
    std::vector<double> axis_r_cos, axis_r_sin, axis_z_cos, axis_z_sin;

    // Unmodified source fields on the first radial_order half surfaces.
    // This is input for a consumer's axis limit, not a Boozer axis ring.
    // Layout [half_surface][source_zeta][source_theta], theta fastest;
    // field ordering and physical-angle units match EquilibriumSnapshot.
    std::vector<double> axis_half_s;
    std::array<std::vector<double>, EquilibriumSnapshot::HALF_FIELD_COUNT>
        axis_half_fields;

    // Present in v6: components in the Boozer basis, sampled on grid.
    // grid.nu and grid.source_theta_shift retain the untruncated angle map.
    std::optional<BoozerMagneticData> magnetic;
};

struct BoozerSurfaceProfiles {
    // Pa, sampled on the same non-axis surfaces as BoozerResult::s/iota.
    std::vector<double> pressure;

    // Covariant coefficients in T m, with physical toroidal angle zeta_b:
    // B = beta*grad(s) + I*grad(theta_b) + G*grad(zeta_b).
    // These are the native angular zero modes of B_theta and B_zeta used
    // by the shift solver, without extra nfp, mu0 or 2*pi factors.
    std::vector<double> current_i;
    std::vector<double> current_g;
};

struct BoozerResult {
    int source_format_version = 0;
    bool lasym = false;
    int source_ns = 0;
    int source_ntheta = 0;
    int source_nzeta = 0;
    int source_mpol = 0;
    int source_ntor = 0;
    int nfp = 0;
    RadialInterpolationOrder radial_order =
        RadialInterpolationOrder::FOUR_POINT;
    double resonance_tolerance = 1.0e-12;

    // Source radial coordinate s=j/(source_ns-1) and iota on non-axis
    // surfaces. s is proportional to toroidal flux only for constant aphi.
    std::vector<double> s;
    std::vector<double> iota;
    BoozerMixedGrid grid;
    MixedGridSpectrum spectrum;

    // Present in new transforms and v5 files. Absent in legacy v3/v4 files:
    // their magnetic-field polarity and axis cannot be recovered from |B|.
    std::optional<BoozerTracingData> tracing;

    // Present in new transforms and v7 files; absent in legacy v3-v6 files.
    std::optional<BoozerSurfaceProfiles> profiles;
};

}  // namespace equilibrium
