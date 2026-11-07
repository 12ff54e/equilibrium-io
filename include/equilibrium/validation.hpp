#pragma once

#include "equilibrium/equilibrium.hpp"
#include "equilibrium/record_validation.hpp"

#include <limits>

namespace equilibrium {

// Structural validation permits absent groups. Consumers choose the groups
// they need; physical validity and convergence remain solver responsibilities.
inline Status validate_group(const EquilibriumView& view, DataGroup group) {
    if (group == DataGroup::GEOMETRY) {
        if (view.ns < 2 || view.mnmax < 1 ||
            (view.families.size() != 6 && view.families.size() != 12))
            return Status("invalid spectral dimensions");
        const auto count = checked_mul(static_cast<std::size_t>(view.ns),
                                       static_cast<std::size_t>(view.mnmax));
        if (!count) return Status("spectral dimensions overflow");
        for (const auto& family : view.families)
            if (family.size() != *count ||
                !std::ranges::all_of(family,
                                     [](double x) { return std::isfinite(x); }))
                return Status("incomplete spectral family");
    }
    if (group == DataGroup::FIELDS && view.fields) {
        const auto& fields = *view.fields;
        if (fields.ntheta < 1 || fields.nzeta < 1)
            return Status("invalid scientific field dimensions");
        const auto plane = checked_mul(static_cast<std::size_t>(fields.ntheta),
                                       static_cast<std::size_t>(fields.nzeta));
        const auto full = plane ? checked_mul(*plane, view.ns) : std::nullopt;
        const auto half =
            plane ? checked_mul(*plane, view.ns - 1) : std::nullopt;
        if (!full || !half) return Status("field dimensions overflow");
        for (const auto& field : fields.half)
            if (field.size() != *half ||
                !std::ranges::all_of(field,
                                     [](double x) { return std::isfinite(x); }))
                return Status("partial half-grid fields");
        for (const auto& field : fields.full)
            if (field.size() != *full ||
                !std::ranges::all_of(field,
                                     [](double x) { return std::isfinite(x); }))
                return Status("partial full-grid fields");
    }
    if (group == DataGroup::PROFILES && view.profiles &&
        !view.profiles->has_half_grid_profiles(view.ns))
        return Status("partial solved radial profiles");
    if (group == DataGroup::GRID && view.grid) {
        const auto& grid = *view.grid;
        if (grid.mpol < 1 || grid.ntor < 0 || grid.nfp < 1 || grid.ntheta < 2 ||
            grid.nzeta < 1 ||
            grid.ntor > std::numeric_limits<int>::max() / grid.nfp)
            return Status("invalid grid metadata");
        const auto modes = checked_mul(static_cast<std::size_t>(grid.mpol),
                                       static_cast<std::size_t>(grid.ntor) + 1);
        if (!modes || *modes != static_cast<std::size_t>(view.mnmax))
            return Status("grid and spectral dimensions differ");
    }
    if (group == DataGroup::INPUT && view.input) {
        const auto& input = view.input->get();
        auto status = validate_input_record(input);
        if (!status) return status;
    }
    if (group == DataGroup::PROVENANCE && view.provenance) {
        const auto& report = view.provenance->get();
        if (static_cast<unsigned>(report.status) > 4 ||
            report.total_effective_iterations < 0 ||
            (report.build.scalar_type != "double" &&
             report.build.scalar_type != "float"))
            return Status("invalid run provenance");
        for (const auto& stage : report.stages) {
            if (stage.ns < 1 || stage.effective_iterations < 0)
                return Status("invalid stage provenance");
            for (const auto& restart : stage.restarts)
                if (restart.iteration < 0)
                    return Status("negative restart iteration");
        }
    }
    return Status();
}

inline Status validate(const EquilibriumView& view) {
    for (auto group :
         {DataGroup::GEOMETRY, DataGroup::FIELDS, DataGroup::GRID,
          DataGroup::PROFILES, DataGroup::INPUT, DataGroup::PROVENANCE}) {
        auto status = validate_group(view, group);
        if (!status) return status;
    }
    if (view.grid && view.fields &&
        (view.fields->ntheta != view.grid->ntheta ||
         view.fields->nzeta != view.grid->nzeta))
        return Status("grid and field dimensions differ");
    if (view.input && view.grid &&
        (*view.grid != grid_metadata(view.input->get()) ||
         view.lasym() != view.input->get().lasym))
        return Status("saved input and scientific metadata differ");
    if (view.input && !view.grid)
        return Status("saved input requires grid metadata");
    return {};
}

}  // namespace equilibrium
