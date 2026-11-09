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

// REAL retains even inactive sin(0) coefficients. Complex reconstruction
// ignores those terms, so every storage mode represents the same geometry.
enum class SpectrumStorage { REAL, COMPLEX, BOTH };

struct BoozerReadOptions {
    // Coordinates are mandatory. Mark consumer prerequisites REQUIRED;
    // optional invalid quantities are omitted with a diagnostic.
    ReadRequest request = ReadRequest::all();
    SpectrumStorage spectrum = SpectrumStorage::REAL;
};

// Writers emit v8; readers accept v7/v8. Complete output requires every
// stored record. Both versions use uniform theta_b and alpha=nfp*zeta, with
// alpha_b=alpha+nfp*nu; nu is in physical toroidal radians.
void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);

BoozerFile read_boozer_binary(const std::filesystem::path& path,
                              const BoozerReadOptions& options = {});

}  // namespace equilibrium

#endif  // EQUILIBRIUM_BOOZER_BINARY_HPP_
