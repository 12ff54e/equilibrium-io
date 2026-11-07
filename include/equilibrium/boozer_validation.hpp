#pragma once

#include "equilibrium/boozer.hpp"
#include "equilibrium/result.hpp"

namespace equilibrium {

// Shared by file readers and consumers of caller-built data. Each quantity
// is checked independently; absent quantities are errors only when requested.
Status validate_boozer_coordinates(const BoozerView& view);
Status validate_quantity(const BoozerView& view, Quantity quantity);

enum class BoozerConsistency {
    FLUX_IOTA,
    SURFACE_CURRENTS,
    MAGNETIC,
    STRENGTH
};
std::span<const Quantity> consistency_quantities(BoozerConsistency check);
// Checks relations only when every participating quantity is present.
// STRENGTH is the legacy B^2*sqrt_g_b identity, not the magnetic-block J.
Status validate_consistency(const BoozerView& view, BoozerConsistency check);

}  // namespace equilibrium
