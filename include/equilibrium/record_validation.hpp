#pragma once

#include "equilibrium/boundary.hpp"
#include "equilibrium/checked_size.hpp"
#include "equilibrium/input_params.hpp"
#include "equilibrium/result.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <string_view>

namespace equilibrium {

// Call after checking equal lengths. Raw ordering and duplicates are retained.
inline std::vector<BoundaryHarmonic> boundary_harmonics(
    std::span<const int> m,
    std::span<const int> n,
    std::span<const double> values) {
    std::vector<BoundaryHarmonic> result;
    result.reserve(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
        result.push_back({m[i], n[i], values[i]});
    return result;
}

// Validate saved scientific structure independently of solver precision,
// initialization, convergence thresholds, and optional vacuum backends.
inline Status validate_input_record(const InputParams& p) {
    if (p.schema != "cumes-config-v1" || p.mpol < 1 || p.ntor < 0 ||
        p.nfp < 1 || p.ntheta < 2 || p.nzeta < 1 ||
        (p.ncurr != 0 && p.ncurr != 1))
        return Status("invalid saved input metadata");
    const auto modes = checked_mul(static_cast<std::size_t>(p.mpol),
                                   static_cast<std::size_t>(p.ntor) + 1);
    if (!modes ||
        *modes > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        p.ntor == std::numeric_limits<int>::max())
        return Status("saved boundary dimensions overflow");
    const auto finite = [](const auto& values) {
        return std::all_of(values.begin(), values.end(),
                           [](double value) { return std::isfinite(value); });
    };
    for (double value : {p.delt, p.phiedge, p.pres_scale, p.adiabatic_index,
                         p.spres_ped, p.bloat, p.curtor, p.tcon0})
        if (!std::isfinite(value))
            return Status("nonfinite saved input scalar");
    if (p.nvacskip < 1) return Status("invalid saved vacuum update interval");
    if (p.embedded_makegrid_parameters) {
        const auto& grid = *p.embedded_makegrid_parameters;
        if (grid.number_of_field_periods != p.nfp ||
            grid.number_of_r_grid_points < 2 ||
            grid.number_of_z_grid_points < 2 ||
            grid.number_of_phi_grid_points < 1 ||
            !std::isfinite(grid.r_grid_minimum) ||
            !std::isfinite(grid.r_grid_maximum) ||
            !std::isfinite(grid.z_grid_minimum) ||
            !std::isfinite(grid.z_grid_maximum) ||
            grid.r_grid_minimum >= grid.r_grid_maximum ||
            grid.z_grid_minimum >= grid.z_grid_maximum)
            return Status("invalid saved MAKEGRID parameters");
    }
    const auto harmonics = [&](const auto& m, const auto& n, const auto& v) {
        if (m.size() != n.size() || m.size() != v.size() || !finite(v))
            return false;
        for (std::size_t i = 0; i < m.size(); ++i)
            if (m[i] < 0 || m[i] >= p.mpol || n[i] < -p.ntor || n[i] > p.ntor)
                return false;
        return true;
    };
    if (!harmonics(p.rbc_m, p.rbc_n, p.rbc_value) ||
        !harmonics(p.zbs_m, p.zbs_n, p.zbs_value) ||
        !harmonics(p.rbs_m, p.rbs_n, p.rbs_value) ||
        !harmonics(p.zbc_m, p.zbc_n, p.zbc_value))
        return Status("invalid saved boundary harmonics");
    for (const auto* values : {&p.rbcc, &p.rbss, &p.zbsc, &p.zbcs})
        if (values->size() != *modes || !finite(*values))
            return Status("incomplete folded boundary");
    for (const auto* values : {&p.rbsc, &p.rbcs, &p.zbcc, &p.zbss})
        if (values->size() != (p.lasym ? *modes : 0) || !finite(*values))
            return Status("inconsistent complementary folded boundary");
    if (!p.lasym && (!p.rbs_value.empty() || !p.zbc_value.empty()))
        return Status("symmetric input contains complementary harmonics");
    const auto axis_size = static_cast<std::size_t>(p.ntor) + 1;
    if (p.raxis_c.size() != axis_size || p.zaxis_s.size() != axis_size ||
        p.raxis_s.size() != (p.lasym ? axis_size : 0) ||
        p.zaxis_c.size() != (p.lasym ? axis_size : 0))
        return Status("incomplete saved magnetic axis");
    for (const auto* values : {&p.raxis_c, &p.zaxis_s, &p.raxis_s, &p.zaxis_c,
                               &p.am, &p.ai, &p.ac, &p.aphi, &p.extcur})
        if (!finite(*values)) return Status("nonfinite saved input array");
    const auto profile = [](std::string_view type, const auto& values,
                            bool allow_two_power) {
        return type == "power_series" ||
               (type == "list" && values.size() >= 2) ||
               (allow_two_power && type == "two_power" && values.size() >= 3);
    };
    if (!profile(p.pmass_type, p.am, true) ||
        !profile(p.piota_type, p.ai, false) ||
        !profile(p.pcurr_type, p.ac, true) ||
        !profile(p.pflux_type, p.aphi, false))
        return Status("invalid saved profile prescription");
    if (p.pflux_type == "list" && std::any_of(p.aphi.begin(), p.aphi.end(),
                                              [](double x) { return x <= 0; }))
        return Status("invalid saved flux-derivative list");
    if (p.pmass_type == "list" &&
        std::any_of(p.am.begin(), p.am.end(), [](double x) { return x < 0; }))
        return Status("invalid saved pressure list");
    int previous_ns = 0;
    if (p.stages.empty()) return Status("missing saved stage schedule");
    for (const auto& stage : p.stages) {
        if (stage.ns < 2 || stage.ns <= previous_ns || stage.max_iter < 1 ||
            !std::isfinite(stage.ftol) || stage.ftol <= 0)
            return Status("invalid saved stage schedule");
        previous_ns = stage.ns;
    }
    struct Boundary {
        int ntor;
        bool lasym;
        std::vector<BoundaryHarmonic> rbc, zbs, rbs, zbc;
    };
    const Boundary raw{p.ntor,
                       p.lasym,
                       boundary_harmonics(p.rbc_m, p.rbc_n, p.rbc_value),
                       boundary_harmonics(p.zbs_m, p.zbs_n, p.zbs_value),
                       boundary_harmonics(p.rbs_m, p.rbs_n, p.rbs_value),
                       boundary_harmonics(p.zbc_m, p.zbc_n, p.zbc_value)};
    const auto folded = fold_boundary(raw, *modes);
    if (folded.rbcc != p.rbcc || folded.rbss != p.rbss ||
        folded.zbsc != p.zbsc || folded.zbcs != p.zbcs ||
        folded.rbsc != p.rbsc || folded.rbcs != p.rbcs ||
        folded.zbcc != p.zbcc || folded.zbss != p.zbss)
        return Status("saved raw and folded input do not match");
    return Status();
}

}  // namespace equilibrium
