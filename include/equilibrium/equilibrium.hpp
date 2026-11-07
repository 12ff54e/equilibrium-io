#pragma once

#include "equilibrium/data_views.hpp"
#include "equilibrium/equilibrium_profiles.hpp"
#include "equilibrium/equilibrium_snapshot.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/read_request.hpp"
#include "equilibrium/run_report.hpp"

#include <functional>
#include <optional>
#include <utility>

namespace equilibrium {

struct GridMetadata {
    int mpol = 0;
    int ntor = 0;
    int nfp = 0;
    int ntheta = 0;
    int nzeta = 0;
    bool operator==(const GridMetadata&) const = default;
};

struct SourceInfo {
    std::optional<OutputFormat> format;
    int version = 0;  // Zero denotes an in-memory producer.
    bool checkpoint = false;
};

template <class Buffer>
struct BasicFields {
    int ntheta = 0;
    int nzeta = 0;
    std::array<Buffer, EquilibriumSnapshot::HALF_FIELD_COUNT> half;
    std::array<Buffer, EquilibriumSnapshot::FULL_FIELD_COUNT> full;
};

// The same scientific layout supports ownership and read-only borrowing.
// Metadata records may likewise be values or reference_wrappers.
template <class Buffer, class Input, class Provenance>
struct BasicEquilibrium {
    int ns = 0;
    int mnmax = 0;
    std::vector<Buffer> families;
    std::optional<BasicFields<Buffer>> fields;
    std::optional<BasicProfiles<Buffer>> profiles;
    std::optional<GridMetadata> grid;
    std::optional<Input> input;
    std::optional<Provenance> provenance;
    SourceInfo source;
    ReadDiagnostics read_diagnostics;
    // Preserve source symmetry when selected reads omit spectral arrays.
    std::optional<bool> symmetry;

    bool has_quantity(Quantity quantity) const {
        const auto nonempty = [](const auto& a) { return !a.empty(); };
        if (quantity == Quantity::GEOMETRY ||
            quantity == Quantity::SPECTRAL_SHIFT) {
            if (families.size() != 6 && families.size() != 12) return false;
            if (lasym() && families.size() != 12) return false;
            for (std::size_t i = 0; i < families.size(); ++i)
                if ((i % 3 == 2) == (quantity == Quantity::SPECTRAL_SHIFT) &&
                    (i < 6 || lasym()) && !nonempty(families[i]))
                    return false;
            return true;
        }
        if (quantity == Quantity::GRID) return grid.has_value();
        if (quantity == Quantity::INPUT) return input.has_value();
        if (quantity == Quantity::PROVENANCE) return provenance.has_value();
        if (!fields) return false;
        if (quantity == Quantity::MAGNETIC_JACOBIAN)
            return nonempty(fields->half[0]);
        if (quantity == Quantity::CONTRAVARIANT_FIELD ||
            quantity == Quantity::COVARIANT_FIELD) {
            const int first = quantity == Quantity::CONTRAVARIANT_FIELD ? 1 : 4;
            for (int i = first; i < first + 3; ++i)
                if (!nonempty(fields->half[i])) return false;
            return true;
        }
        if (quantity == Quantity::CURRENT_DENSITY) {
            for (const auto& values : fields->full)
                if (!nonempty(values)) return false;
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
    bool has_fields() const { return has_group(DataGroup::FIELDS); }
    bool has_grid() const { return grid.has_value(); }
    bool has_input() const { return input.has_value(); }
    bool has_provenance() const { return provenance.has_value(); }
    bool has_profiles() const { return profiles.has_value(); }

    bool lasym() const { return symmetry.value_or(families.size() == 12); }

    int field_periods() const { return grid.value().nfp; }
    SpectralGrid spectral_grid() const {
        return {ns, grid.value().mpol, grid.value().ntor,
                SpectralLayout::MODE_MAJOR};
    }
    double full_s(std::size_t radial) const {
        return full_radial_s(ns, radial);
    }
    double half_s(std::size_t radial) const {
        return half_radial_s(ns, radial);
    }

    std::optional<VectorFieldView> contravariant_field() const& {
        return vector_field(Quantity::CONTRAVARIANT_FIELD, 1);
    }
    std::optional<VectorFieldView> covariant_field() const& {
        return vector_field(Quantity::COVARIANT_FIELD, 4);
    }
    std::optional<ScalarFieldView> magnetic_jacobian() const& {
        if (!has_quantity(Quantity::MAGNETIC_JACOBIAN)) return std::nullopt;
        return ScalarFieldView{{ns - 1, fields->ntheta, fields->nzeta},
                               fields->half[0]};
    }
    std::optional<CurrentDensityView> current_density() const& {
        if (!has_quantity(Quantity::CURRENT_DENSITY)) return std::nullopt;
        const FieldGrid shape{ns, fields->ntheta, fields->nzeta};
        return CurrentDensityView{
            {shape, {fields->full[0], fields->full[1], fields->full[2]}},
            {shape, {fields->full[3], fields->full[4], fields->full[5]}}};
    }
    std::optional<GeometryCoefficientsView> geometry_coefficients() const& {
        if (!grid || !has_geometry()) return std::nullopt;
        const auto family = [&](std::size_t i) -> std::span<const double> {
            return i < families.size() ? std::span<const double>(families[i])
                                       : std::span<const double>{};
        };
        return GeometryCoefficientsView{
            spectral_grid(), lasym(),   family(0), family(3), family(6),
            family(9),       family(1), family(4), family(7), family(10)};
    }
    std::optional<ShiftCoefficientsView> lambda_coefficients() const& {
        if (!grid || !has_quantity(Quantity::SPECTRAL_SHIFT))
            return std::nullopt;
        return ShiftCoefficientsView{
            spectral_grid(),
            lasym(),
            families[2],
            families[5],
            lasym() ? std::span<const double>(families[8])
                    : std::span<const double>{},
            lasym() ? std::span<const double>(families[11])
                    : std::span<const double>{}};
    }
    std::optional<AxisGeometryView> axis_geometry() const& {
        const auto geometry = geometry_coefficients();
        if (!geometry) return std::nullopt;
        const auto axis = [&](std::span<const double> values) {
            return StridedValues{values,
                                 static_cast<std::size_t>(grid->ntor + 1),
                                 static_cast<std::size_t>(ns)};
        };
        return AxisGeometryView{grid->nfp,
                                {axis(geometry->rcc), axis(geometry->rcs),
                                 axis(geometry->zcc), axis(geometry->zcs)}};
    }

    std::optional<VectorFieldView> contravariant_field() const&& = delete;
    std::optional<VectorFieldView> covariant_field() const&& = delete;
    std::optional<ScalarFieldView> magnetic_jacobian() const&& = delete;
    std::optional<CurrentDensityView> current_density() const&& = delete;
    std::optional<GeometryCoefficientsView> geometry_coefficients() const&& =
        delete;
    std::optional<ShiftCoefficientsView> lambda_coefficients() const&& = delete;
    std::optional<AxisGeometryView> axis_geometry() const&& = delete;

   private:
    std::optional<VectorFieldView> vector_field(Quantity quantity,
                                                int first) const {
        if (!has_quantity(quantity)) return std::nullopt;
        return VectorFieldView{{ns - 1, fields->ntheta, fields->nzeta},
                               {fields->half[first], fields->half[first + 1],
                                fields->half[first + 2]}};
    }
};

using Equilibrium =
    BasicEquilibrium<std::vector<double>, InputParams, RunProvenance>;
using EquilibriumView =
    BasicEquilibrium<std::span<const double>,
                     std::reference_wrapper<const InputParams>,
                     std::reference_wrapper<const RunProvenance>>;

inline GridMetadata grid_metadata(const InputParams& input) {
    return {input.mpol, input.ntor, input.nfp, input.ntheta, input.nzeta};
}

// Owners must outlive the returned view and keep their buffers intact.
inline EquilibriumView as_view(const Equilibrium& owner) {
    EquilibriumView view;
    view.ns = owner.ns;
    view.mnmax = owner.mnmax;
    view.families.assign(owner.families.begin(), owner.families.end());
    view.grid = owner.grid;
    view.source = owner.source;
    view.read_diagnostics = owner.read_diagnostics;
    view.symmetry = owner.symmetry;
    if (owner.fields) {
        auto& fields = view.fields.emplace();
        fields.ntheta = owner.fields->ntheta;
        fields.nzeta = owner.fields->nzeta;
        for (std::size_t i = 0; i < fields.half.size(); ++i)
            fields.half[i] = owner.fields->half[i];
        for (std::size_t i = 0; i < fields.full.size(); ++i)
            fields.full[i] = owner.fields->full[i];
    }
    if (owner.profiles) view.profiles = as_view(*owner.profiles);
    if (owner.input) view.input = std::cref(*owner.input);
    if (owner.provenance) view.provenance = std::cref(*owner.provenance);
    return view;
}
EquilibriumView as_view(Equilibrium&&) = delete;
EquilibriumView as_view(const Equilibrium&&) = delete;

// Compose independent host owners. This adapter does not validate or invent
// missing metadata. A partial field block remains present so validation fails.
inline EquilibriumView as_equilibrium_view(const EquilibriumSnapshot& state) {
    EquilibriumView view;
    view.ns = state.ns;
    view.mnmax = state.mnmax;
    view.families.assign(state.families.begin(), state.families.end());
    view.symmetry = state.lasym();
    view.source.version = state.source_version;
    view.source.checkpoint = state.source_checkpoint;
    if (state.has_any_derived_fields()) {
        auto& fields = view.fields.emplace();
        fields.ntheta = state.ntheta;
        fields.nzeta = state.nzeta;
        for (std::size_t i = 0; i < fields.half.size(); ++i)
            fields.half[i] = state.half_fields[i];
        for (std::size_t i = 0; i < fields.full.size(); ++i)
            fields.full[i] = state.full_fields[i];
    }
    return view;
}
EquilibriumView as_equilibrium_view(EquilibriumSnapshot&&) = delete;
EquilibriumView as_equilibrium_view(const EquilibriumSnapshot&&) = delete;

inline SnapshotView as_snapshot_view(const EquilibriumView& view) {
    SnapshotView state;
    state.ns = view.ns;
    state.mnmax = view.mnmax;
    state.families = view.families;
    if (!view.lasym() && state.families.size() == 12) state.families.resize(6);
    state.source_version = view.source.version;
    state.source_checkpoint = view.source.checkpoint;
    if (view.fields) {
        state.ntheta = view.fields->ntheta;
        state.nzeta = view.fields->nzeta;
        state.half_fields = view.fields->half;
        state.full_fields = view.fields->full;
    }
    return state;
}

// Transfer decoded buffers without copying their numerical contents.
inline Equilibrium take_equilibrium(EquilibriumSnapshot state) {
    Equilibrium owner;
    owner.ns = state.ns;
    owner.mnmax = state.mnmax;
    owner.symmetry = state.lasym();
    owner.families = std::move(state.families);
    owner.source.version = state.source_version;
    owner.source.checkpoint = state.source_checkpoint;
    if (state.has_any_derived_fields()) {
        owner.fields = {state.ntheta, state.nzeta, std::move(state.half_fields),
                        std::move(state.full_fields)};
    }
    return owner;
}

}  // namespace equilibrium
