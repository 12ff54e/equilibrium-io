#pragma once

#include "equilibrium/result.hpp"
#include "equilibrium/spdata.hpp"

#include <string>

namespace equilibrium {

// Write prepared coefficients in MCT's spdata text format. The path may have
// any suffix. Values retain the caller's units and are rounded to E18.10.
// Validate before opening output; failures before rename preserve the target.
// A directory-sync error after rename reports unconfirmed durability.
// Link equilibrium::io; no optional codec or interpolation library is needed.
Status write_spdata(const std::string& path, const SpdataView& data);
inline Status write_spdata(const std::string& path, const Spdata& data) {
    return write_spdata(path, as_view(data));
}

}  // namespace equilibrium
