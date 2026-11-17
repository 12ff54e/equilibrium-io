#pragma once

#include "equilibrium/checked_size.hpp"
#include "equilibrium/result.hpp"
#include "equilibrium/spdata.hpp"

#include <array>
#include <cmath>

namespace equilibrium {

// Check the coefficient record, not force balance or spline accuracy.
inline Status validate_spdata(const SpdataView& data) {
    // MCT writes I4 dimensions but reads them as whitespace-separated values.
    // A four-digit lst would join the preceding lsp token.
    if (data.lsp < 2 || data.lsp > 9999)
        return std::string("spdata lsp must be in [2, 9999]");
    if (data.lst < 1 || data.lst > 999)
        return std::string("spdata lst must be in [1, 999]");
    if (!std::isfinite(data.psi_wall) || data.psi_wall <= 0.0 ||
        !std::isfinite(data.psi_sep) || data.psi_sep <= 0.0)
        return std::string("spdata flux limits must be finite and positive");
    if (data.meta_data.size() > (1u << 20) ||
        data.meta_data.find_first_of(std::string_view("\r\n\0", 3)) !=
            std::string_view::npos)
        return std::string("spdata metadata must be one line of at most 1 MiB");

    const auto radial_coefficients = checked_mul(data.lsp, 9);
    const auto field_size = radial_coefficients
                                ? checked_mul(*radial_coefficients, data.lst)
                                : std::nullopt;
    const auto profile_size = checked_mul(data.lsp, 3);
    if (!field_size || !profile_size)
        return std::string("spdata coefficient count overflow");
    const std::array fields{data.b, data.r, data.z, data.j};
    const std::array profiles{data.q, data.f,       data.i,
                              data.p, data.r_minor, data.psi_t};
    constexpr std::array field_names{"b", "r", "z", "j"};
    constexpr std::array profile_names{"q", "f", "i", "p", "r_minor", "psi_t"};
    const auto check_values = [](std::span<const double> values,
                                 std::size_t expected,
                                 const char* name) -> Status {
        if (values.size() != expected)
            return "spdata " + std::string(name) + " has the wrong size";
        for (double value : values)
            if (!std::isfinite(value))
                return "spdata " + std::string(name) + " must be finite";
        return {};
    };
    for (std::size_t k = 0; k < fields.size(); ++k) {
        if (auto status = check_values(fields[k], *field_size, field_names[k]);
            !status)
            return status;
        // The axis is one point: its value cannot depend on theta.
        for (std::size_t t = 0; t < data.lst; ++t) {
            if (fields[k][t] != fields[k][0] ||
                fields[k][data.field_index(0, 3, t)] != 0.0 ||
                fields[k][data.field_index(0, 6, t)] != 0.0)
                return "spdata " + std::string(field_names[k]) +
                       " has an angle-dependent axis value";
        }
    }
    for (std::size_t k = 0; k < profiles.size(); ++k)
        if (auto status =
                check_values(profiles[k], *profile_size, profile_names[k]);
            !status)
            return status;
    if (data.r.front() <= 0.0)
        return std::string("spdata axis radius must be positive");
    return {};
}

inline Status validate_spdata(const Spdata& data) {
    return validate_spdata(as_view(data));
}

}  // namespace equilibrium
