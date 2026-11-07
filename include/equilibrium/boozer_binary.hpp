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
    // Core coordinates/dimensions are always read. All other groups follow
    // this request. The default requires geometry/strength and requests the
    // optional tracing, magnetic/maps and surface-profile groups.
    ReadRequest groups = [] {
        ReadRequest request;
        request[DataGroup::GEOMETRY] = Requirement::REQUIRED;
        request[DataGroup::STRENGTH] = Requirement::REQUIRED;
        request[Quantity::ROTATIONAL_TRANSFORM] = Requirement::REQUIRED;
        for (auto group :
             {DataGroup::TRACING, DataGroup::MAGNETIC, DataGroup::PROFILES})
            request[group] = Requirement::OPTIONAL;
        return request;
    }();
    SpectrumStorage spectrum = SpectrumStorage::COMPLEX;
};

// New transforms write v7 (fields/maps, pressure, I/G), or v8 for flux
// derivative lists. Legacy results retain their v3-v6 format when rewritten.
// All versions use uniform theta_b and unchanged alpha=nfp*zeta, with
// alpha_b=alpha+nfp*nu; nu is in physical toroidal radians.
void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);

BoozerFile read_boozer_binary(const std::filesystem::path& path);
BoozerFile read_boozer_binary(const std::filesystem::path& path,
                              const BoozerReadOptions& options);

}  // namespace equilibrium

#endif  // EQUILIBRIUM_BOOZER_BINARY_HPP_
