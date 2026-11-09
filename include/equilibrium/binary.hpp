#pragma once

#include "equilibrium/equilibrium.hpp"

namespace equilibrium {

enum class ProfileWritePolicy { REQUIRE_PERSISTENCE, OMIT };

struct WriteOptions {
    ProfileWritePolicy profiles = ProfileWritePolicy::REQUIRE_PERSISTENCE;
};

// Native results and checkpoints share the quantity-request API.
Result<Equilibrium> read_binary(
    const std::string& path,
    const ReadRequest& request = ReadRequest::all());

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
