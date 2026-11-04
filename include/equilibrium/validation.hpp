#pragma once

#include "equilibrium/equilibrium.hpp"
#include "equilibrium/record_validation.hpp"

#include <limits>

namespace equilibrium {

// Structural validation permits absent groups. Consumers choose the groups
// they need; physical validity and convergence remain solver responsibilities.
inline Status validate(const EquilibriumView& view) {
    if (view.ns < 2 || view.mnmax < 1 ||
        (view.families.size() != 6 && view.families.size() != 12))
        return Status("invalid spectral dimensions");
    const auto count = checked_mul(static_cast<std::size_t>(view.ns),
                                   static_cast<std::size_t>(view.mnmax));
    if (!count) return Status("spectral dimensions overflow");
    for (const auto& family : view.families)
        if (family.size() != *count)
            return Status("incomplete spectral family");
    if (view.fields) {
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
            if (field.size() != *half)
                return Status("partial half-grid fields");
        for (const auto& field : fields.full)
            if (field.size() != *full)
                return Status("partial full-grid fields");
    }
    if (view.profiles && !view.profiles->has_half_grid_profiles(view.ns))
        return Status("partial solved radial profiles");
    if (view.grid) {
        const auto& grid = *view.grid;
        if (grid.mpol < 1 || grid.ntor < 0 || grid.nfp < 1 || grid.ntheta < 2 ||
            grid.nzeta < 1 ||
            grid.ntor > std::numeric_limits<int>::max() / grid.nfp)
            return Status("invalid grid metadata");
        const auto modes = checked_mul(static_cast<std::size_t>(grid.mpol),
                                       static_cast<std::size_t>(grid.ntor) + 1);
        if (!modes || *modes != static_cast<std::size_t>(view.mnmax))
            return Status("grid and spectral dimensions differ");
        if (view.fields && (view.fields->ntheta != grid.ntheta ||
                            view.fields->nzeta != grid.nzeta))
            return Status("grid and field dimensions differ");
    }
    if (view.input) {
        const auto& input = view.input->get();
        auto status = validate_input_record(input);
        if (!status) return status;
        if (!view.grid || *view.grid != grid_metadata(input) ||
            view.lasym() != input.lasym)
            return Status("saved input and scientific metadata differ");
    }
    if (view.provenance) {
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

}  // namespace equilibrium
