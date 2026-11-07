#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <string>

namespace equilibrium {

enum class Requirement { REQUIRED, OPTIONAL, IGNORED };
enum class Availability { AVAILABLE, ABSENT, NOT_REQUESTED, INVALID, PARTIAL };
enum class Quantity {
    GEOMETRY,
    SPECTRAL_SHIFT,
    GRID,
    CONTRAVARIANT_FIELD,
    COVARIANT_FIELD,
    MAGNETIC_JACOBIAN,
    CURRENT_DENSITY,
    INPUT,
    PROVENANCE,
    FIELD_STRENGTH,
    STRENGTH_JACOBIAN,
    B2_JACOBIAN,
    ROTATIONAL_TRANSFORM,
    FLUX_INPUT,
    FLUX_DERIVATIVES,
    AXIS_GEOMETRY,
    AXIS_FIELDS,
    TOROIDAL_SHIFT,
    SOURCE_THETA_SHIFT,
    PRESSURE,
    SURFACE_CURRENTS,
    COUNT
};
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

inline constexpr Quantity half_field_quantity(std::size_t component) {
    return component == 0  ? Quantity::MAGNETIC_JACOBIAN
           : component < 4 ? Quantity::CONTRAVARIANT_FIELD
                           : Quantity::COVARIANT_FIELD;
}
inline constexpr Quantity native_family_quantity(std::size_t family) {
    return family % 3 == 2 ? Quantity::SPECTRAL_SHIFT : Quantity::GEOMETRY;
}

inline std::span<const Quantity> group_quantities(DataGroup group) {
    using Q = Quantity;
    static constexpr std::array geometry{Q::GEOMETRY, Q::SPECTRAL_SHIFT};
    static constexpr std::array fields{Q::MAGNETIC_JACOBIAN,
                                       Q::CONTRAVARIANT_FIELD,
                                       Q::COVARIANT_FIELD, Q::CURRENT_DENSITY};
    static constexpr std::array grid{Q::GRID};
    static constexpr std::array input{Q::INPUT};
    static constexpr std::array provenance{Q::PROVENANCE};
    static constexpr std::array profiles{Q::PRESSURE, Q::SURFACE_CURRENTS};
    static constexpr std::array strength{Q::FIELD_STRENGTH,
                                         Q::STRENGTH_JACOBIAN, Q::B2_JACOBIAN};
    static constexpr std::array tracing{Q::FLUX_INPUT, Q::FLUX_DERIVATIVES,
                                        Q::AXIS_GEOMETRY, Q::AXIS_FIELDS};
    static constexpr std::array magnetic{
        Q::MAGNETIC_JACOBIAN, Q::CONTRAVARIANT_FIELD, Q::COVARIANT_FIELD,
        Q::TOROIDAL_SHIFT, Q::SOURCE_THETA_SHIFT};
    switch (group) {
        case DataGroup::GEOMETRY:
            return geometry;
        case DataGroup::FIELDS:
            return fields;
        case DataGroup::GRID:
            return grid;
        case DataGroup::INPUT:
            return input;
        case DataGroup::PROVENANCE:
            return provenance;
        case DataGroup::PROFILES:
            return profiles;
        case DataGroup::STRENGTH:
            return strength;
        case DataGroup::TRACING:
            return tracing;
        case DataGroup::MAGNETIC:
            return magnetic;
        default:
            return {};
    }
}

// Groups are assignment presets, not atomic validation units. Assignments
// expand immediately; a later quantity or group assignment wins. Optional
// quantities fail independently. Framing needed to locate required data must
// still be intact even when an intervening quantity is ignored.
struct ReadRequest {
    std::array<Requirement, static_cast<std::size_t>(Quantity::COUNT)>
        quantities;

    ReadRequest() { quantities.fill(Requirement::IGNORED); }
    Requirement& operator[](Quantity quantity) {
        return quantities[static_cast<std::size_t>(quantity)];
    }
    Requirement operator[](Quantity quantity) const {
        return quantities[static_cast<std::size_t>(quantity)];
    }
    struct GroupPreset {
        ReadRequest& request;
        DataGroup group;
        GroupPreset& operator=(Requirement requirement) {
            for (auto q : group_quantities(group)) request[q] = requirement;
            return *this;
        }
    };
    GroupPreset operator[](DataGroup group) { return {*this, group}; }
    // The strongest requirement among the preset's members.
    Requirement strongest_requirement(DataGroup group) const {
        auto result = Requirement::IGNORED;
        for (auto q : group_quantities(group))
            if ((*this)[q] < result) result = (*this)[q];
        return result;
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
    std::array<GroupStatus, static_cast<std::size_t>(Quantity::COUNT)>
        quantities;
    GroupStatus& operator[](Quantity q) {
        return quantities[static_cast<std::size_t>(q)];
    }
    const GroupStatus& operator[](Quantity q) const {
        return quantities[static_cast<std::size_t>(q)];
    }
    GroupStatus& operator[](DataGroup group) {
        return groups[static_cast<std::size_t>(group)];
    }
    const GroupStatus& operator[](DataGroup group) const {
        return groups[static_cast<std::size_t>(group)];
    }
};

inline const char* quantity_name(Quantity quantity) {
    if (static_cast<std::size_t>(quantity) >=
        static_cast<std::size_t>(Quantity::COUNT))
        return "unsupported quantity";
    constexpr std::array names{"geometry",
                               "spectral shift",
                               "grid",
                               "contravariant field",
                               "covariant field",
                               "magnetic Jacobian",
                               "current density",
                               "input",
                               "provenance",
                               "field strength",
                               "strength Jacobian",
                               "B2 Jacobian",
                               "rotational transform",
                               "flux input",
                               "flux derivatives",
                               "axis geometry",
                               "axis fields",
                               "toroidal shift",
                               "source theta shift",
                               "pressure",
                               "surface currents"};
    return names[static_cast<std::size_t>(quantity)];
}

inline const char* group_name(DataGroup group) {
    if (static_cast<std::size_t>(group) >=
        static_cast<std::size_t>(DataGroup::COUNT))
        return "unsupported group";
    constexpr std::array names{"geometry", "fields",     "grid",
                               "input",    "provenance", "profiles",
                               "strength", "tracing",    "magnetic"};
    return names[static_cast<std::size_t>(group)];
}

}  // namespace equilibrium
