#ifndef EQUILIBRIUM_BOOZER_OUTPUT_HPP_
#define EQUILIBRIUM_BOOZER_OUTPUT_HPP_

#include "equilibrium/boozer.hpp"

#include <cstdint>
#include <filesystem>
#include <string_view>

namespace equilibrium {

enum class BoozerOutputFormat : std::uint8_t {
    BINARY = 0,
    NETCDF = 1,
    HDF5 = 2,
};

struct BoozerOutputSpec {
    BoozerOutputFormat format = BoozerOutputFormat::BINARY;
    std::filesystem::path path;
};

// Resolve .bin, .nc, .h5, or .hdf5 case-insensitively. Unknown and missing
// suffixes are rejected instead of receiving misleading binary contents.
BoozerOutputSpec resolve_boozer_output_spec(const std::filesystem::path& path);

bool boozer_output_format_available(BoozerOutputFormat format);
std::string_view boozer_output_suffix(BoozerOutputFormat format);

// Write the selected container, using v8 when flux lists are present.
// Publish atomically; failures before rename preserve the destination.
// The result can also remain in memory without serialization.
void write_boozer_output(const BoozerOutputSpec& spec,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);

inline void write_boozer_output(const std::filesystem::path& path,
                                const BoozerResult& result,
                                const std::filesystem::path& source_path) {
    write_boozer_output(resolve_boozer_output_spec(path), result, source_path);
}

}  // namespace equilibrium

#endif  // EQUILIBRIUM_BOOZER_OUTPUT_HPP_
