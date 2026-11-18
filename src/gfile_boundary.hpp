#pragma once

#include "equilibrium/gfile.hpp"

namespace equilibrium::gfile_detail {
// Profile-only input. The caller supplies a fresh boundary before writing it.
GfileInput read_profiles(const std::filesystem::path& path,
                         const GfileInputSettings& settings);
// Input must have no boundary harmonics. Return maximum sampled error in m.
double fit_boundary(InputParams& input,
                    std::span<const LcfsPoint> sampled,
                    LcfsPoint axis,
                    double minor_radius,
                    double tolerance);
}  // namespace equilibrium::gfile_detail
