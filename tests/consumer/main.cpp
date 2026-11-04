#include <equilibrium/equilibrium.hpp>

#ifdef EQUILIBRIUM_COMPONENT_binary
#include <equilibrium/binary.hpp>
#endif
#ifdef EQUILIBRIUM_COMPONENT_io
#include <equilibrium/io.hpp>
#endif

int main() {
    equilibrium::Equilibrium owner;
    const auto view = equilibrium::as_view(owner);
    if (!view.families.empty()) return 1;
#ifdef EQUILIBRIUM_COMPONENT_binary
    if (equilibrium::read_binary("does-not-exist.bin")) return 2;
#endif
#ifdef EQUILIBRIUM_COMPONENT_io
    if (equilibrium::read("does-not-exist.bin")) return 2;
#endif
}
