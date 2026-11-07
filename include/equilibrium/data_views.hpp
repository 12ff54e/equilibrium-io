#pragma once

#include "equilibrium/checked_size.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <span>

namespace equilibrium {

// Borrowed views never copy numerical arrays. Their owner must remain alive
// and must not resize or replace the buffers while a view is in use.
struct FieldGrid {
    int radial_points = 0;
    int ntheta = 0;
    int nzeta = 0;

    std::optional<std::size_t> size() const {
        if (radial_points < 1 || ntheta < 1 || nzeta < 1) return std::nullopt;
        auto plane = checked_mul(static_cast<std::size_t>(ntheta),
                                 static_cast<std::size_t>(nzeta));
        return plane ? checked_mul(*plane, radial_points) : std::nullopt;
    }
    // Theta is contiguous. Indices must lie within the grid.
    std::size_t index(std::size_t radial, int theta, int zeta) const {
        return (radial * nzeta + zeta) * ntheta + theta;
    }
};

struct VectorFieldView {
    FieldGrid grid;
    std::array<std::span<const double>, 3> components;
};

struct ScalarFieldView {
    FieldGrid grid;
    std::span<const double> values;
};

struct CurrentDensityView {
    VectorFieldView contravariant;
    VectorFieldView covariant;
};

struct FluxDerivativesView {
    // Signed flux/(2*pi) derivatives, Wb/rad per unit source s.
    std::span<const double> toroidal;
    std::span<const double> poloidal;
    double s(std::size_t radial) const {
        return static_cast<double>(radial) / (toroidal.size() - 1);
    }
};

struct SurfaceCurrentsView {
    std::span<const double> i;  // T m, physical toroidal angle.
    std::span<const double> g;
};

struct ScaledProfileView {
    std::span<const double> values;
    double divisor = 1;
    std::size_t size() const { return values.size(); }
    double operator[](std::size_t index) const {
        return values[index] / divisor;
    }
};

inline constexpr double BOOZER_ANGULAR_INTEGRAL =
    (2 * std::numbers::pi) * (2 * std::numbers::pi);

// The legacy angular mean B^2 J, not the signed Jacobian in magnetic data.
inline ScaledProfileView boozer_b2_jacobian(std::span<const double> b2j00) {
    return {b2j00, BOOZER_ANGULAR_INTEGRAL};
}

inline double full_radial_s(int source_ns, std::size_t radial) {
    return static_cast<double>(radial) / (source_ns - 1);
}
inline double half_radial_s(int source_ns, std::size_t radial) {
    return (radial + .5) / (source_ns - 1);
}

enum class SpectralLayout { MODE_MAJOR, SURFACE_MAJOR };
struct SpectralGrid {
    int surfaces = 0;
    int mpol = 0;
    int ntor = 0;
    SpectralLayout layout = SpectralLayout::MODE_MAJOR;

    // Native: m-major modes, surface contiguous. Boozer: n-major modes,
    // mode contiguous. Indices must lie within the spectral grid.
    std::size_t index(std::size_t surface, int m, int n) const {
        if (layout == SpectralLayout::MODE_MAJOR)
            return static_cast<std::size_t>(m * (ntor + 1) + n) * surfaces +
                   surface;
        return (surface * (ntor + 1) + n) * mpol + m;
    }
};

struct GeometryCoefficientsView {
    SpectralGrid grid;
    bool lasym = false;
    std::span<const double> rcc, rss, rsc, rcs;
    std::span<const double> zsc, zcs, zcc, zss;
};

struct ShiftCoefficientsView {
    SpectralGrid grid;
    bool lasym = false;
    std::span<const double> sc, cs, cc, ss;
};

struct StridedValues {
    std::span<const double> values;
    std::size_t count = 0;
    std::size_t stride = 1;
    double operator[](std::size_t index) const {
        return values.empty() ? 0.0 : values[index * stride];
    }
};

struct AxisGeometryView {
    int nfp = 1;
    std::array<StridedValues, 4> families;  // Rcos, Rsin, Zcos, Zsin.
    std::array<double, 2> cylindrical(double phi) const {
        std::array<double, 2> result{};
        for (std::size_t n = 0; n < families[0].count; ++n) {
            const double phase = static_cast<double>(n) * nfp * phi;
            const double c = std::cos(phase), s = std::sin(phase);
            result[0] += families[0][n] * c + families[1][n] * s;
            result[1] += families[2][n] * c + families[3][n] * s;
        }
        return result;
    }
};

}  // namespace equilibrium
