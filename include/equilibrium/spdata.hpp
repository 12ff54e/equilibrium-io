#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace equilibrium {

// Prepared quadratic coefficients for an axisymmetric equilibrium. Coordinate
// conversion and spline fitting belong to the producer. See docs/spdata.md.
template <class Buffer, class String>
struct BasicSpdata {
    String meta_data;
    std::size_t lsp = 0;
    std::size_t lst = 0;
    double psi_wall = 0.0;
    double psi_sep = 0.0;

    // [radial record][coefficient][theta], with lst unique angles.
    // coefficient = radial_order + 3 * theta_order; both orders run 0..2.
    // The first radial record uses sqrt(psi); later records use delta_psi.
    Buffer b, r, z, j;

    // [radial record][coefficient], three coefficients per record.
    // Only r_minor uses sqrt(psi) in the first record; the others use psi.
    Buffer q, f, i, p, r_minor, psi_t;

    std::size_t field_index(std::size_t radial,
                            std::size_t coefficient,
                            std::size_t theta) const {
        return (radial * 9 + coefficient) * lst + theta;
    }
    std::size_t profile_index(std::size_t radial,
                              std::size_t coefficient) const {
        return radial * 3 + coefficient;
    }
};

using Spdata = BasicSpdata<std::vector<double>, std::string>;
using SpdataView = BasicSpdata<std::span<const double>, std::string_view>;

inline SpdataView as_view(const Spdata& owner) {
    return {owner.meta_data, owner.lsp,     owner.lst,  owner.psi_wall,
            owner.psi_sep,   owner.b,       owner.r,    owner.z,
            owner.j,         owner.q,       owner.f,    owner.i,
            owner.p,         owner.r_minor, owner.psi_t};
}
SpdataView as_view(Spdata&&) = delete;
SpdataView as_view(const Spdata&&) = delete;

}  // namespace equilibrium
