// equilibrium_profiles.hpp — host radial profiles associated with a solved
// equilibrium. These are in-memory solver outputs, separate from optimizer
// observables and target/reduction policy.
#ifndef EQUILIBRIUM_INCLUDE_IO_EQUILIBRIUM_PROFILES_HPP_
#define EQUILIBRIUM_INCLUDE_IO_EQUILIBRIUM_PROFILES_HPP_

#include <cstddef>
#include <span>
#include <vector>

namespace equilibrium {

template <class Buffer>
struct BasicProfiles {
    // Half-grid profiles at s=(j+1/2)/(ns-1). Flux derivatives use the public
    // VMEC wout convention in webers per unit s: the internal oriented
    // phip/chip values are multiplied by sign(J)*2*pi during capture.
    Buffer toroidal_flux_derivative;
    Buffer poloidal_flux_derivative;
    Buffer rotational_transform;

    // VMEC-compatible covariant flux functions on the half grid:
    //   I(s) = <B_theta>  (wout buco)
    //   G(s) = <B_zeta>   (wout bvco).
    // The angle brackets are the normalized VMEC angular quadrature, not a
    // physical-area-weighted average.
    Buffer poloidal_covariant_field;
    Buffer toroidal_covariant_field;

    bool has_half_grid_profiles(int ns) const {
        if (ns < 2) return false;
        const std::size_t expected = static_cast<std::size_t>(ns - 1);
        return toroidal_flux_derivative.size() == expected &&
               poloidal_flux_derivative.size() == expected &&
               rotational_transform.size() == expected &&
               poloidal_covariant_field.size() == expected &&
               toroidal_covariant_field.size() == expected;
    }
};

using EquilibriumProfiles = BasicProfiles<std::vector<double>>;
using ProfilesView = BasicProfiles<std::span<const double>>;
inline ProfilesView as_view(const EquilibriumProfiles& p) {
    return {p.toroidal_flux_derivative, p.poloidal_flux_derivative,
            p.rotational_transform, p.poloidal_covariant_field,
            p.toroidal_covariant_field};
}
ProfilesView as_view(EquilibriumProfiles&&) = delete;
ProfilesView as_view(const EquilibriumProfiles&&) = delete;

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INCLUDE_IO_EQUILIBRIUM_PROFILES_HPP_
