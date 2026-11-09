#pragma once

#include "equilibrium/equilibrium_snapshot.hpp"
#include "equilibrium/result.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace equilibrium {

struct SpectralPayload {
    int ns = 0;
    int mnmax = 0;
    int components = 6;
    std::vector<std::uint8_t> bytes;
};

// Inspect only the spectral portion of cuMES result binaries (versions 10/11).
// These functions accept files without derived fields or provenance. They do
// not validate the rest of the file or the finiteness of the coefficients.
// Checkpoints and scientific containers use read_binary/read instead.
// The bytes API preserves signed zero and NaN payloads for bitwise comparison.
Result<SpectralPayload> read_spectral_payload(const std::string& path);
Result<EquilibriumSnapshot> read_spectral_state(const std::string& path);

}  // namespace equilibrium
