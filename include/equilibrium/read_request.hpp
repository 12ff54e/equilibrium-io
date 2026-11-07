#pragma once

#include <array>
#include <cstddef>
#include <string>

namespace equilibrium {

enum class Requirement { REQUIRED, OPTIONAL, IGNORED };
enum class Availability { AVAILABLE, ABSENT, NOT_REQUESTED, INVALID };
enum class DataGroup {
    GEOMETRY,
    FIELDS,
    GRID,
    INPUT,
    PROVENANCE,
    PROFILES,
    STRENGTH,
    TRACING,
    MAGNETIC,
    COUNT
};

// A request starts with every group ignored. Required groups must be usable;
// optional groups report missing or invalid data without failing the read.
// Structural damage that prevents locating requested data still fails a read.
struct ReadRequest {
    std::array<Requirement, static_cast<std::size_t>(DataGroup::COUNT)> groups;

    ReadRequest() { groups.fill(Requirement::IGNORED); }
    Requirement& operator[](DataGroup group) {
        return groups[static_cast<std::size_t>(group)];
    }
    Requirement operator[](DataGroup group) const {
        return groups[static_cast<std::size_t>(group)];
    }
};

struct GroupStatus {
    Availability availability = Availability::NOT_REQUESTED;
    std::string diagnostic;
};

// These describe the read, even if the caller later changes the public data.
// has_* queries on an equilibrium describe its current exposed groups.
struct ReadDiagnostics {
    std::array<GroupStatus, static_cast<std::size_t>(DataGroup::COUNT)> groups;
    GroupStatus& operator[](DataGroup group) {
        return groups[static_cast<std::size_t>(group)];
    }
    const GroupStatus& operator[](DataGroup group) const {
        return groups[static_cast<std::size_t>(group)];
    }
};

inline const char* group_name(DataGroup group) {
    constexpr std::array names{"geometry", "fields",     "grid",
                               "input",    "provenance", "profiles",
                               "strength", "tracing",    "magnetic"};
    return names[static_cast<std::size_t>(group)];
}

}  // namespace equilibrium
