#pragma once

#include "equilibrium/data_views.hpp"
#include "equilibrium/equilibrium_snapshot.hpp"
#include "equilibrium/read_request.hpp"

#include <array>
#include <complex>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace equilibrium {

// Boozer records use the existing mixed grid: theta_b is uniform and
// alpha=nfp*zeta remains the source field-period angle. These are data records;
// coordinate transforms and interpolation operations belong to consumers.
enum class RadialInterpolationOrder {
    TWO_POINT = 2,
    FOUR_POINT = 4,
};

template <class Buffer>
struct BasicBoozerMixedGrid {
    int source_ns = 0;
    int first_surface = 1;
    int ntheta = 0;
    int nzeta = 0;

    Buffer source_theta_p;
    Buffer source_theta;
    // Smooth periodic displacement theta_native - theta_b, without wrapping.
    // Together with nu, this specifies the full source-to-Boozer map.
    Buffer source_theta_shift;
    Buffer r;
    Buffer z;
    Buffer nu;
    Buffer b;
    Buffer sqrtg_b;
    Buffer b2j00;
};
using BoozerMixedGrid = BasicBoozerMixedGrid<std::vector<double>>;

// Original real parity coefficients: [surface][mode], n-major, m >= 0, n >= 0.
// REAL reads preserve these values exactly, including inactive zero modes.
template <class Buffer, class Indices>
struct BasicRealBoozerSpectrum {
    bool lasym = false;
    int mmax = 0;
    int nmax = 0;
    Indices m;
    Indices n;
    Buffer rmncc;
    Buffer rmnss;
    Buffer zmnsc;
    Buffer zmncs;
    Buffer numnsc;
    Buffer numncs;
    Buffer rmnsc;
    Buffer rmncs;
    Buffer zmncc;
    Buffer zmnss;
    Buffer numncc;
    Buffer numnss;
};
using RealBoozerSpectrum =
    BasicRealBoozerSpectrum<std::vector<double>, std::vector<int>>;

template <class ComplexBuffer, class Indices>
struct BasicMixedGridSpectrum {
    int source_ns = 0;
    int first_surface = 1;
    int mmax = 0;
    int nmax = 0;

    // Mode ordering is n-major, then m from -mmax through mmax.
    Indices m;
    Indices n;
    ComplexBuffer r;
    ComplexBuffer z;
    ComplexBuffer nu;
};
using MixedGridSpectrum =
    BasicMixedGridSpectrum<std::vector<std::complex<double>>, std::vector<int>>;

template <class Buffer>
struct BasicBoozerMagneticData {
    // Signed Jacobian and six components in (s, theta_b, physical zeta_b).
    // Ordering matches EquilibriumSnapshot::HALF_FIELD_COUNT; sampling matches
    // BoozerMixedGrid: [non-axis surface][source alpha][theta_b].
    std::array<Buffer, EquilibriumSnapshot::HALF_FIELD_COUNT> fields;
};
using BoozerMagneticData = BasicBoozerMagneticData<std::vector<double>>;

template <class Buffer, class String>
struct BasicBoozerTracingData {
    // Original cuMES flux-profile inputs (SI), without a polarity change.
    double phiedge = 0.0;
    Buffer aphi;
    String pflux_type = "power_series";

    // Signed <sqrt(g) B^zeta> and <sqrt(g) B^theta>, interpolated to all
    // source full-grid surfaces s=j/(source_ns-1), including the axis.
    // zeta is the physical toroidal angle; these are flux/(2*pi) derivatives.
    Buffer toroidal_flux_derivative;
    Buffer poloidal_flux_derivative;

    // Converged m=0 axis geometry, all n=0..source_ntor, independent of the
    // output truncation: R/Z=sum[cos[n]*cos(n*alpha)+sin[n]*sin(n*alpha)].
    Buffer axis_r_cos, axis_r_sin, axis_z_cos, axis_z_sin;

    // Unmodified source fields on the first radial_order half surfaces.
    // This is input for a consumer's axis limit, not a Boozer axis ring.
    // Layout [half_surface][source_zeta][source_theta], theta fastest;
    // field ordering and physical-angle units match EquilibriumSnapshot.
    Buffer axis_half_s;
    std::array<Buffer, EquilibriumSnapshot::HALF_FIELD_COUNT> axis_half_fields;

    // Present in v6: components in the Boozer basis, sampled on grid.
    // grid.nu and grid.source_theta_shift retain the untruncated angle map.
    std::optional<BasicBoozerMagneticData<Buffer>> magnetic;
};
using BoozerTracingData =
    BasicBoozerTracingData<std::vector<double>, std::string>;

template <class Buffer>
struct BasicBoozerSurfaceProfiles {
    // Pa, sampled on the same non-axis surfaces as BoozerResult::s/iota.
    Buffer pressure;

    // Covariant coefficients in T m, with physical toroidal angle zeta_b:
    // B = beta*grad(s) + I*grad(theta_b) + G*grad(zeta_b).
    // These are the native angular zero modes of B_theta and B_zeta used
    // by the shift solver, without extra nfp, mu0 or 2*pi factors.
    Buffer current_i;
    Buffer current_g;
};
using BoozerSurfaceProfiles = BasicBoozerSurfaceProfiles<std::vector<double>>;

template <class Buffer, class Indices, class ComplexBuffer, class String>
struct BasicBoozerResult {
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
    Buffer s;
    Buffer iota;
    BasicBoozerMixedGrid<Buffer> grid;
    BasicMixedGridSpectrum<ComplexBuffer, Indices> spectrum;
    std::optional<BasicRealBoozerSpectrum<Buffer, Indices>> real_spectrum;
    ReadDiagnostics read_diagnostics;

    // Present in new transforms and v5 files. Absent in legacy v3/v4 files:
    // their magnetic-field polarity and axis cannot be recovered from |B|.
    std::optional<BasicBoozerTracingData<Buffer, String>> tracing;

    // Present in new transforms and v7 files; absent in legacy v3-v6 files.
    std::optional<BasicBoozerSurfaceProfiles<Buffer>> profiles;

    // For selected reads with tracing ignored, magnetic data remain usable.
    // Complete/legacy reads retain the existing tracing->magnetic location.
    std::optional<BasicBoozerMagneticData<Buffer>> independent_magnetic;

    bool has_quantity(Quantity q) const {
        if (q == Quantity::GEOMETRY || q == Quantity::SPECTRAL_SHIFT) {
            if (real_spectrum) {
                const auto& r = *real_spectrum;
                if (q == Quantity::SPECTRAL_SHIFT)
                    return !r.numnsc.empty() && !r.numncs.empty() &&
                           (!lasym || (!r.numncc.empty() && !r.numnss.empty()));
                return !r.rmncc.empty() && !r.rmnss.empty() &&
                       !r.zmnsc.empty() && !r.zmncs.empty() &&
                       (!lasym || (!r.rmnsc.empty() && !r.rmncs.empty() &&
                                   !r.zmncc.empty() && !r.zmnss.empty()));
            }
            return q == Quantity::GEOMETRY
                       ? !spectrum.r.empty() && !spectrum.z.empty()
                       : !spectrum.nu.empty();
        }
        if (q == Quantity::GRID)
            return source_ns >= 2 && grid.ntheta > 0 && grid.nzeta > 0;
        if (q == Quantity::ROTATIONAL_TRANSFORM) return !iota.empty();
        if (q == Quantity::FIELD_STRENGTH) return !grid.b.empty();
        if (q == Quantity::STRENGTH_JACOBIAN) return !grid.sqrtg_b.empty();
        if (q == Quantity::B2_JACOBIAN) return !grid.b2j00.empty();
        if (q == Quantity::TOROIDAL_SHIFT) return !grid.nu.empty();
        if (q == Quantity::SOURCE_THETA_SHIFT)
            return !grid.source_theta_shift.empty();
        if (q == Quantity::PRESSURE)
            return profiles && !profiles->pressure.empty();
        if (q == Quantity::SURFACE_CURRENTS)
            return profiles && !profiles->current_i.empty() &&
                   !profiles->current_g.empty();
        if (q == Quantity::FLUX_INPUT) return tracing && !tracing->aphi.empty();
        if (q == Quantity::FLUX_DERIVATIVES)
            return tracing && !tracing->toroidal_flux_derivative.empty() &&
                   !tracing->poloidal_flux_derivative.empty();
        if (q == Quantity::AXIS_GEOMETRY)
            return tracing && !tracing->axis_r_cos.empty() &&
                   !tracing->axis_r_sin.empty() &&
                   !tracing->axis_z_cos.empty() && !tracing->axis_z_sin.empty();
        if (q == Quantity::AXIS_FIELDS) {
            if (!tracing || tracing->axis_half_s.empty()) return false;
            for (const auto& values : tracing->axis_half_fields)
                if (values.empty()) return false;
            return true;
        }
        const auto data = magnetic_data();
        if (!data) return false;
        if (q == Quantity::MAGNETIC_JACOBIAN)
            return !data->get().fields[0].empty();
        if (q == Quantity::CONTRAVARIANT_FIELD ||
            q == Quantity::COVARIANT_FIELD) {
            const int first = q == Quantity::CONTRAVARIANT_FIELD ? 1 : 4;
            for (int i = first; i < first + 3; ++i)
                if (data->get().fields[i].empty()) return false;
            return true;
        }
        return false;
    }
    bool has_group(DataGroup group) const {
        for (auto q : group_quantities(group))
            if (!has_quantity(q)) return false;
        return group != DataGroup::COUNT;
    }
    bool has_geometry() const { return has_quantity(Quantity::GEOMETRY); }
    bool has_real_spectrum() const { return real_spectrum.has_value(); }
    bool has_complex_spectrum() const {
        return !spectrum.r.empty() || !spectrum.nu.empty();
    }
    bool has_strength() const { return has_group(DataGroup::STRENGTH); }
    bool has_tracing() const { return has_group(DataGroup::TRACING); }
    bool has_profiles() const { return has_group(DataGroup::PROFILES); }
    std::optional<std::reference_wrapper<const BasicBoozerMagneticData<Buffer>>>
    magnetic_data() const& {
        if (independent_magnetic) return std::cref(*independent_magnetic);
        if (tracing && tracing->magnetic) return std::cref(*tracing->magnetic);
        return std::nullopt;
    }
    bool has_magnetic() const { return has_group(DataGroup::MAGNETIC); }

    int field_periods() const { return nfp; }
    double full_s(std::size_t r) const { return full_radial_s(source_ns, r); }
    double half_s(std::size_t r) const { return half_radial_s(source_ns, r); }
    SpectralGrid spectral_grid() const {
        return {source_ns - grid.first_surface,
                real_spectrum ? real_spectrum->mmax + 1 : spectrum.mmax + 1,
                real_spectrum ? real_spectrum->nmax : spectrum.nmax,
                SpectralLayout::SURFACE_MAJOR};
    }
    std::optional<GeometryCoefficientsView> geometry_coefficients() const& {
        if (!real_spectrum || !has_geometry()) return std::nullopt;
        const auto& r = *real_spectrum;
        return GeometryCoefficientsView{
            spectral_grid(), lasym,   r.rmncc, r.rmnss, r.rmnsc,
            r.rmncs,         r.zmnsc, r.zmncs, r.zmncc, r.zmnss};
    }
    std::optional<ShiftCoefficientsView> nu_coefficients() const& {
        if (!real_spectrum || !has_quantity(Quantity::SPECTRAL_SHIFT))
            return std::nullopt;
        const auto& r = *real_spectrum;
        return ShiftCoefficientsView{spectral_grid(), lasym,    r.numnsc,
                                     r.numncs,        r.numncc, r.numnss};
    }
    std::optional<VectorFieldView> contravariant_field() const& {
        return vector_field(Quantity::CONTRAVARIANT_FIELD, 1);
    }
    std::optional<VectorFieldView> covariant_field() const& {
        return vector_field(Quantity::COVARIANT_FIELD, 4);
    }
    std::optional<ScalarFieldView> magnetic_jacobian() const& {
        if (!has_quantity(Quantity::MAGNETIC_JACOBIAN)) return std::nullopt;
        return ScalarFieldView{field_grid(), magnetic_data()->get().fields[0]};
    }
    std::optional<ScalarFieldView> field_strength() const& {
        if (!has_quantity(Quantity::FIELD_STRENGTH)) return std::nullopt;
        return ScalarFieldView{field_grid(), grid.b};
    }
    std::optional<ScalarFieldView> strength_jacobian() const& {
        if (!has_quantity(Quantity::STRENGTH_JACOBIAN)) return std::nullopt;
        return ScalarFieldView{field_grid(), grid.sqrtg_b};
    }
    std::optional<ScaledProfileView> b2_jacobian() const& {
        if (!has_quantity(Quantity::B2_JACOBIAN)) return std::nullopt;
        return boozer_b2_jacobian(grid.b2j00);
    }
    std::span<const double> rotational_transform() const& { return iota; }
    std::span<const double> toroidal_shift() const& { return grid.nu; }
    std::span<const double> source_theta_shift() const& {
        return grid.source_theta_shift;
    }
    std::span<const double> pressure() const& {
        return profiles ? std::span<const double>(profiles->pressure)
                        : std::span<const double>{};
    }
    std::optional<FluxDerivativesView> flux_derivatives() const& {
        if (!has_quantity(Quantity::FLUX_DERIVATIVES)) return std::nullopt;
        return FluxDerivativesView{tracing->toroidal_flux_derivative,
                                   tracing->poloidal_flux_derivative};
    }
    std::optional<SurfaceCurrentsView> surface_currents() const& {
        if (!has_quantity(Quantity::SURFACE_CURRENTS)) return std::nullopt;
        return SurfaceCurrentsView{profiles->current_i, profiles->current_g};
    }
    std::optional<AxisGeometryView> axis_geometry() const& {
        if (!has_quantity(Quantity::AXIS_GEOMETRY)) return std::nullopt;
        const auto axis = [](std::span<const double> v) {
            return StridedValues{v, v.size(), 1};
        };
        return AxisGeometryView{
            nfp,
            {axis(tracing->axis_r_cos), axis(tracing->axis_r_sin),
             axis(tracing->axis_z_cos), axis(tracing->axis_z_sin)}};
    }
    std::optional<std::reference_wrapper<const BasicBoozerMagneticData<Buffer>>>
    magnetic_data() const&& = delete;
    std::optional<GeometryCoefficientsView> geometry_coefficients() const&& =
        delete;
    std::optional<ShiftCoefficientsView> nu_coefficients() const&& = delete;
    std::optional<VectorFieldView> contravariant_field() const&& = delete;
    std::optional<VectorFieldView> covariant_field() const&& = delete;
    std::optional<ScalarFieldView> magnetic_jacobian() const&& = delete;
    std::optional<ScalarFieldView> field_strength() const&& = delete;
    std::optional<ScalarFieldView> strength_jacobian() const&& = delete;
    std::optional<ScaledProfileView> b2_jacobian() const&& = delete;
    std::span<const double> rotational_transform() const&& = delete;
    std::span<const double> toroidal_shift() const&& = delete;
    std::span<const double> source_theta_shift() const&& = delete;
    std::span<const double> pressure() const&& = delete;
    std::optional<FluxDerivativesView> flux_derivatives() const&& = delete;
    std::optional<SurfaceCurrentsView> surface_currents() const&& = delete;
    std::optional<AxisGeometryView> axis_geometry() const&& = delete;

   private:
    FieldGrid field_grid() const {
        return {source_ns - grid.first_surface, grid.ntheta, grid.nzeta};
    }
    std::optional<VectorFieldView> vector_field(Quantity q, int first) const {
        if (!has_quantity(q)) return std::nullopt;
        const auto& f = magnetic_data()->get().fields;
        return VectorFieldView{field_grid(),
                               {f[first], f[first + 1], f[first + 2]}};
    }
};
using BoozerResult = BasicBoozerResult<std::vector<double>,
                                       std::vector<int>,
                                       std::vector<std::complex<double>>,
                                       std::string>;

using BoozerView = BasicBoozerResult<std::span<const double>,
                                     std::span<const int>,
                                     std::span<const std::complex<double>>,
                                     std::string_view>;

inline BoozerView as_view(const BoozerResult& owner) {
    BoozerView view;
    view.read_diagnostics = owner.read_diagnostics;
    view.source_format_version = owner.source_format_version;
    view.lasym = owner.lasym;
    view.source_ns = owner.source_ns;
    view.source_ntheta = owner.source_ntheta;
    view.source_nzeta = owner.source_nzeta;
    view.source_mpol = owner.source_mpol;
    view.source_ntor = owner.source_ntor;
    view.nfp = owner.nfp;
    view.radial_order = owner.radial_order;
    view.resonance_tolerance = owner.resonance_tolerance;
    view.s = owner.s;
    view.iota = owner.iota;
    view.grid = {owner.grid.source_ns,
                 owner.grid.first_surface,
                 owner.grid.ntheta,
                 owner.grid.nzeta,
                 owner.grid.source_theta_p,
                 owner.grid.source_theta,
                 owner.grid.source_theta_shift,
                 owner.grid.r,
                 owner.grid.z,
                 owner.grid.nu,
                 owner.grid.b,
                 owner.grid.sqrtg_b,
                 owner.grid.b2j00};
    const auto& spectrum = owner.spectrum;
    view.spectrum = {spectrum.source_ns, spectrum.first_surface,
                     spectrum.mmax,      spectrum.nmax,
                     spectrum.m,         spectrum.n,
                     spectrum.r,         spectrum.z,
                     spectrum.nu};
    if (owner.real_spectrum) {
        const auto& r = *owner.real_spectrum;
        view.real_spectrum = {r.lasym,  r.mmax,  r.nmax,  r.m,     r.n,
                              r.rmncc,  r.rmnss, r.zmnsc, r.zmncs, r.numnsc,
                              r.numncs, r.rmnsc, r.rmncs, r.zmncc, r.zmnss,
                              r.numncc, r.numnss};
    }
    if (owner.tracing) {
        const auto& t = *owner.tracing;
        auto& v = view.tracing.emplace();
        v.phiedge = t.phiedge;
        v.aphi = t.aphi;
        v.pflux_type = t.pflux_type;
        v.toroidal_flux_derivative = t.toroidal_flux_derivative;
        v.poloidal_flux_derivative = t.poloidal_flux_derivative;
        v.axis_r_cos = t.axis_r_cos;
        v.axis_r_sin = t.axis_r_sin;
        v.axis_z_cos = t.axis_z_cos;
        v.axis_z_sin = t.axis_z_sin;
        v.axis_half_s = t.axis_half_s;
        for (std::size_t i = 0; i < v.axis_half_fields.size(); ++i)
            v.axis_half_fields[i] = t.axis_half_fields[i];
        if (t.magnetic) {
            auto& m = v.magnetic.emplace();
            for (std::size_t i = 0; i < m.fields.size(); ++i)
                m.fields[i] = t.magnetic->fields[i];
        }
    }
    if (owner.independent_magnetic) {
        auto& m = view.independent_magnetic.emplace();
        for (std::size_t i = 0; i < m.fields.size(); ++i)
            m.fields[i] = owner.independent_magnetic->fields[i];
    }
    if (owner.profiles)
        view.profiles = {owner.profiles->pressure, owner.profiles->current_i,
                         owner.profiles->current_g};
    return view;
}
BoozerView as_view(BoozerResult&&) = delete;
BoozerView as_view(const BoozerResult&&) = delete;

}  // namespace equilibrium
