#ifndef EQUILIBRIUM_BOOZER_BINARY_HPP_
#define EQUILIBRIUM_BOOZER_BINARY_HPP_

#include "equilibrium/boozer.hpp"

#include <filesystem>
#include <string>

namespace equilibrium {

// Check dimensions, finite data and optional-group consistency; throws on
// error.
void validate_boozer_result(const BoozerResult& result);

struct BoozerFile {
    std::string source_path;
    BoozerResult result;
};

// New transforms write v7 (fields/maps, pressure, I/G), or v8 for flux
// derivative lists. Legacy results retain their v3-v6 format when rewritten.
// All versions use uniform theta_b and unchanged alpha=nfp*zeta, with
// alpha_b=alpha+nfp*nu; nu is in physical toroidal radians.
void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);

BoozerFile read_boozer_binary(const std::filesystem::path& path);

}  // namespace equilibrium

#endif  // EQUILIBRIUM_BOOZER_BINARY_HPP_
