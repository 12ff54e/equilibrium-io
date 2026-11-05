#ifndef EQUILIBRIUM_LIST_PROFILE_HPP_
#define EQUILIBRIUM_LIST_PROFILE_HPP_

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string_view>

namespace equilibrium::profile_detail {

// Uniform samples on [0,1]; callers validate the extent and sample values.
inline double evaluate_list(std::span<const double> values,
                            double x,
                            bool integrate = false) {
    if (!std::isfinite(x)) return x;
    if (x <= 0.0) return integrate ? 0.0 : values.front();
    if (x >= 1.0) {
        if (!integrate) return values.back();
        x = 1.0;
    }
    const double intervals = static_cast<double>(values.size() - 1);
    const double position = x * intervals;
    const auto i =
        std::min(static_cast<std::size_t>(position), values.size() - 2);
    const double fraction = position - static_cast<double>(i);
    const double value = std::lerp(values[i], values[i + 1], fraction);
    if (!integrate) return value;
    double integral = 0.0;
    for (std::size_t j = 0; j < i; ++j)
        integral += std::lerp(values[j], values[j + 1], 0.5) / intervals;
    return integral + fraction / intervals * std::lerp(values[i], value, 0.5);
}

inline void validate_flux_type(std::string_view type,
                               std::span<const double> values) {
    if (type == "power_series") return;
    if (type != "list" || values.size() < 2 ||
        !std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value) && value > 0.0;
        }))
        throw std::invalid_argument(
            "flux profile requires power_series or a list of at least two "
            "positive finite derivatives");
}

}  // namespace equilibrium::profile_detail

#endif  // EQUILIBRIUM_LIST_PROFILE_HPP_
