#pragma once

#include <cstddef>
#include <cstdlib>
#include <vector>

namespace equilibrium {
struct BoundaryHarmonic {
    int m = 0;
    int n = 0;
    double value = 0.0;
};

struct FoldedBoundary {
    std::vector<double> rbcc;  // R: cos(mθ)cos(nζ)
    std::vector<double> rbss;  // R: sin(mθ)sin(nζ)
    std::vector<double> zbsc;  // Z: sin(mθ)cos(nζ)
    std::vector<double> zbcs;  // Z: cos(mθ)sin(nζ)
    std::vector<double> rbsc;  // R: sin(mθ)cos(nζ), lasym only
    std::vector<double> rbcs;  // R: cos(mθ)sin(nζ), lasym only
    std::vector<double> zbcc;  // Z: cos(mθ)cos(nζ), lasym only
    std::vector<double> zbss;  // Z: sin(mθ)sin(nζ), lasym only

    std::size_t size() const { return rbcc.size(); }
};

// Fold raw signed-n boundary harmonics into the n>=0 product basis (blueprint
// §4.2), matching the legacy foldBoundary exactly:
//   rbcc[m,n] = rbc[m,+n] + rbc[m,-n]        (all m)
//   rbss[m,n] = rbc[m,+n] - rbc[m,-n]        (m > 0)
//   zbsc[m,n] = zbs[m,+n] + zbs[m,-n]        (m > 0)
//   zbcs[m,n] = zbs[m,-n] - zbs[m,+n]        (all m)
template <class Boundary>
FoldedBoundary fold_boundary(const Boundary& spec, std::size_t modes) {
    FoldedBoundary fb;
    fb.rbcc.assign(modes, 0.0);
    fb.rbss.assign(modes, 0.0);
    fb.zbsc.assign(modes, 0.0);
    fb.zbcs.assign(modes, 0.0);
    const std::size_t ntorp1 = static_cast<std::size_t>(spec.ntor + 1);
    auto idx = [&](int m, int n) {
        return static_cast<std::size_t>(m) * ntorp1 +
               static_cast<std::size_t>(n);
    };
    for (const auto& e : spec.rbc) {
        const int n = std::abs(e.n);
        fb.rbcc[idx(e.m, n)] += e.value;
        const double sn = (e.n > 0) ? 1.0 : (e.n < 0) ? -1.0 : 0.0;
        if (e.m > 0) fb.rbss[idx(e.m, n)] += sn * e.value;
    }
    for (const auto& e : spec.zbs) {
        const int n = std::abs(e.n);
        const double sn = (e.n > 0) ? 1.0 : (e.n < 0) ? -1.0 : 0.0;
        if (e.m > 0) fb.zbsc[idx(e.m, n)] += e.value;
        fb.zbcs[idx(e.m, n)] -= sn * e.value;
    }
    if (spec.lasym) {
        fb.rbsc.assign(modes, 0.0);
        fb.rbcs.assign(modes, 0.0);
        fb.zbcc.assign(modes, 0.0);
        fb.zbss.assign(modes, 0.0);
        for (const auto& e : spec.rbs) {
            const int n = std::abs(e.n);
            const double sn = (e.n > 0) ? 1.0 : (e.n < 0) ? -1.0 : 0.0;
            if (e.m > 0) fb.rbsc[idx(e.m, n)] += e.value;
            fb.rbcs[idx(e.m, n)] -= sn * e.value;
        }
        for (const auto& e : spec.zbc) {
            const int n = std::abs(e.n);
            const double sn = (e.n > 0) ? 1.0 : (e.n < 0) ? -1.0 : 0.0;
            fb.zbcc[idx(e.m, n)] += e.value;
            if (e.m > 0) fb.zbss[idx(e.m, n)] += sn * e.value;
        }
    }
    return fb;
}

}  // namespace equilibrium
