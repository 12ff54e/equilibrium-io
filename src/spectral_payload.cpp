#include "equilibrium/spectral_payload.hpp"

#include "equilibrium/writer_helpers.hpp"
#include "format_versions.hpp"

#include <array>
#include <bit>
#include <fstream>
#include <string_view>

namespace equilibrium {
namespace {

std::int32_t read_i32(const std::array<std::uint8_t, 24>& bytes,
                      std::size_t offset) {
    std::uint32_t value = 0;
    for (std::size_t i = 0; i < 4; ++i)
        value |= std::uint32_t(bytes[offset + i]) << (8 * i);
    return std::bit_cast<std::int32_t>(value);
}

}  // namespace

Result<SpectralPayload> read_spectral_payload(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return Result<SpectralPayload>("could not open " + path);
    const auto end = file.tellg();
    if (end < 8)
        return Result<SpectralPayload>("truncated state header: " + path);
    const auto file_size = static_cast<std::uint64_t>(end);
    file.seekg(0);
    std::array<std::uint8_t, 24> header{};
    file.read(reinterpret_cast<char*>(header.data()), header.size());
    const auto header_size = file.gcount();
    const bool versioned =
        std::string_view(reinterpret_cast<const char*>(header.data()), 8) ==
        "CUMES001";
    SpectralPayload payload;
    std::size_t offset = 8;
    if (versioned) {
        if (header_size < 20)
            return Result<SpectralPayload>("truncated state header: " + path);
        const int version = read_i32(header, 8);
        if (!io_detail::supported_version(version, io_detail::NATIVE_VERSION))
            return Result<SpectralPayload>("unsupported state version: " +
                                           path);
        payload.ns = read_i32(header, 12);
        payload.mnmax = read_i32(header, 16);
        offset = 24;
        if (header_size < static_cast<std::streamsize>(offset))
            return Result<SpectralPayload>("truncated component count: " +
                                           path);
        payload.components = read_i32(header, 20);
    } else {
        return Result<SpectralPayload>("not a cuMES state container: " + path);
    }
    if (payload.ns < 1 || payload.mnmax < 1 ||
        (payload.components != 6 && payload.components != 12))
        return Result<SpectralPayload>("invalid spectral dimensions: " + path);
    const auto family_size =
        checked_mul(static_cast<std::size_t>(payload.ns),
                    static_cast<std::size_t>(payload.mnmax));
    if (!family_size || *family_size > io_detail::MAX_STATE_ELEMENTS_PER_FAMILY)
        return Result<SpectralPayload>(
            "spectral dimensions exceed resource limit: " + path);
    const auto size =
        checked_mul(*family_size, payload.components * sizeof(double));
    if (!size || file_size < offset || *size > file_size - offset)
        return Result<SpectralPayload>("truncated state payload: " + path);
    file.clear();
    file.seekg(static_cast<std::streamoff>(offset));
    payload.bytes.resize(*size);
    file.read(reinterpret_cast<char*>(payload.bytes.data()),
              payload.bytes.size());
    if (!file)
        return Result<SpectralPayload>("could not read state payload: " + path);
    return payload;
}

Result<EquilibriumSnapshot> read_spectral_state(const std::string& path) {
    auto result = read_spectral_payload(path);
    if (!result) return Result<EquilibriumSnapshot>(result.error());
    const auto& payload = result.value();
    EquilibriumSnapshot state;
    state.ns = payload.ns;
    state.mnmax = payload.mnmax;
    state.families.resize(payload.components);
    std::size_t offset = 0;
    for (auto& family : state.families) {
        family.resize(state.family_size());
        for (auto& value : family) {
            std::uint64_t bits = 0;
            for (std::size_t i = 0; i < sizeof(double); ++i)
                bits |= std::uint64_t(payload.bytes[offset++]) << (8 * i);
            value = std::bit_cast<double>(bits);
        }
    }
    return state;
}

}  // namespace equilibrium
