#pragma once

#include "equilibrium/equilibrium.hpp"

namespace equilibrium {

enum class ReadMode { COMPLETE, STATE_ONLY };
enum class ProfileWritePolicy { REQUIRE_PERSISTENCE, OMIT };

struct WriteOptions {
    ProfileWritePolicy profiles = ProfileWritePolicy::REQUIRE_PERSISTENCE;
};

// COMPLETE requires saved input; result files also require their full report.
// STATE_ONLY deliberately skips provenance, including unsupported historical
// trailers, and returns no input, grid metadata, or run provenance.
Result<Equilibrium> read_binary(const std::string& path,
                                ReadMode mode = ReadMode::COMPLETE);
Result<Equilibrium> read_binary(const std::string& path,
                                const ReadRequest& request);

enum class BinaryKind { NATIVE, BOOZER };
// Classification belongs to the codecs; successful detection does not imply
// that the file contains the groups a consumer needs.
Result<BinaryKind> detect_binary_kind(const std::string& path);
Status write_binary(const std::string& path,
                    const EquilibriumView& view,
                    WriteOptions options = {});
inline Status write_binary(const std::string& path,
                           const Equilibrium& owner,
                           WriteOptions options = {}) {
    return write_binary(path, as_view(owner), options);
}

}  // namespace equilibrium
