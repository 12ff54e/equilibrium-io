#include "equilibrium/spectral_payload.hpp"
#include "test_support.hpp"

#include <array>
#include <bit>
#include <fstream>
#include <limits>

using namespace equilibrium;
using test::check;

namespace {
void append(std::vector<std::uint8_t>& data, std::uint64_t value, int bytes) {
    for (int i = 0; i < bytes; ++i) data.push_back((value >> (8 * i)) & 255);
}
void save(const std::vector<std::uint8_t>& data) {
    std::ofstream out("payload.bin", std::ios::binary);
    out.write(reinterpret_cast<const char*>(data.data()), data.size());
}
std::vector<std::uint8_t> fixture(int version, int components = 6) {
    std::vector<std::uint8_t> data;
    if (version) {
        for (char c : std::string_view("CUMES001")) data.push_back(c);
        append(data, version, 4);
    }
    append(data, 2, 4);
    append(data, 1, 4);
    if (version >= 9) append(data, components, 4);
    constexpr std::array<std::uint64_t, 4> bits = {
        0, 0x8000000000000000, 0x7ff8000000001234, 0x3ff0000000000000};
    for (int i = 0; i < 2 * components; ++i) append(data, bits[i % 4], 8);
    return data;
}
}  // namespace

int main() {
    for (int version : {0, 1, 8, 9, 10, 11}) {
        const int components = version >= 9 ? 12 : 6;
        auto bytes = fixture(version, components);
        const std::size_t offset = version == 0 ? 8 : version >= 9 ? 24 : 20;
        const std::vector<std::uint8_t> expected(bytes.begin() + offset,
                                                 bytes.end());
        save(bytes);
        auto payload = read_spectral_payload("payload.bin", version == 0);
        check(payload && payload.value().bytes == expected,
              "spectral bytes survive without derived fields or provenance");
        auto state = read_spectral_state("payload.bin", version == 0);
        check(state && state.value().ns == 2 && state.value().mnmax == 1 &&
                  state.value().components() == components,
              "numeric state uses the same payload layout");
        if (state) {
            check(std::bit_cast<std::uint64_t>(state.value().families[0][1]) ==
                          0x8000000000000000 &&
                      std::bit_cast<std::uint64_t>(
                          state.value().families[1][0]) == 0x7ff8000000001234,
                  "decoding preserves signed zero and NaN payload bits");
        }
        if (!version) {
            check(!read_spectral_payload("payload.bin"),
                  "legacy requires opt-in");
            bytes.push_back(0);
            save(bytes);
            check(!read_spectral_payload("payload.bin", true),
                  "legacy requires exact length");
        } else {
            bytes.push_back(42);
            save(bytes);
            payload = read_spectral_payload("payload.bin");
            check(payload && payload.value().bytes == expected,
                  "trailing bytes do not affect bitwise comparison");
        }
        bytes.resize(offset + expected.size() - 1);
        save(bytes);
        check(!read_spectral_payload("payload.bin", version == 0),
              "short payload rejected");
    }
    for (std::size_t size : {0u, 7u, 19u, 20u, 23u}) {
        auto bytes = fixture(10, 12);
        bytes.resize(size);
        save(bytes);
        check(!read_spectral_payload("payload.bin"), "short header rejected");
    }
    auto bytes = fixture(12, 12);
    save(bytes);
    check(!read_spectral_payload("payload.bin"), "future version rejected");
    bytes = fixture(10, 7);
    save(bytes);
    check(!read_spectral_payload("payload.bin"), "wrong family count rejected");
    bytes = fixture(1);
    bytes[12] = bytes[13] = bytes[14] = 255;
    bytes[15] = 127;
    save(bytes);
    check(!read_spectral_payload("payload.bin"),
          "oversized dimensions rejected before allocation");
    bytes[12] = bytes[13] = bytes[14] = bytes[15] = 0;
    save(bytes);
    check(!read_spectral_payload("payload.bin"), "zero dimensions rejected");
    std::filesystem::remove("payload.bin");
    return test::summary();
}
