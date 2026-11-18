#pragma once

#include "equilibrium/gfile.hpp"

namespace equilibrium {

// Both endpoints are included when reachable by the signed step.
struct MillerRange {
    double first, last, step;
};

struct MillerScan {
    // Unset ranges take a single measured value from the source LCFS.
    std::optional<MillerRange> r0, kappa_upper, kappa_lower, delta_upper,
        delta_lower;
    std::size_t max_cases = 10000;
};

// Write one cuMES input per Cartesian combination and an index.csv with
// requested and fitted LCFS statistics. Requires a new output directory.
// The source supplies minor radius, Z reference and full-edge profiles.
std::size_t write_miller_batch(const std::filesystem::path& gfile,
                               const std::filesystem::path& output_directory,
                               const MillerScan& scan,
                               const GfileInputSettings& settings = {});

}  // namespace equilibrium
