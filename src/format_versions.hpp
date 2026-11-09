#pragma once

namespace equilibrium::io_detail {

inline constexpr int NATIVE_VERSION = 11;
inline constexpr int CHECKPOINT_VERSION = 9;
inline constexpr int BOOZER_VERSION = 8;

constexpr bool supported_version(int version, int current) {
    return version == current || version == current - 1;
}

}  // namespace equilibrium::io_detail
