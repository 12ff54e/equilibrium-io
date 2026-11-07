#pragma once

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

    bool has_geometry() const { return !families.empty(); }
    bool has_fields() const { return fields.has_value(); }
    bool has_grid() const { return grid.has_value(); }
    bool has_input() const { return input.has_value(); }
    bool has_provenance() const { return provenance.has_value(); }
    bool has_profiles() const { return profiles.has_value(); }

    bool lasym() const { return families.size() == 12; }
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
