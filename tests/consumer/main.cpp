#include <equilibrium/equilibrium.hpp>
#include <equilibrium/spdata_validation.hpp>

#ifdef EQUILIBRIUM_COMPONENT_binary
#include <equilibrium/binary.hpp>
#include <equilibrium/boozer_binary.hpp>
#include <equilibrium/spectral_payload.hpp>
#endif
#ifdef EQUILIBRIUM_COMPONENT_io
#include <equilibrium/io.hpp>
#include <equilibrium/spdata_output.hpp>
#endif

#ifdef EQUILIBRIUM_COMPONENT_gfile
#include <sstream>

#include <equilibrium/miller_batch.hpp>
#endif

int main() {
#ifdef EQUILIBRIUM_COMPONENT_gfile
    const std::vector<equilibrium::LcfsPoint> points{
        {4, 0}, {3, 2}, {2, 0}, {3, -1}};
    const auto shape = equilibrium::calculate_lcfs_shape(points, 0);
    if (shape.elongation_upper != 2 || shape.elongation_lower != 1) return 7;
    try {
        equilibrium::write_miller_batch("does-not-exist.g", "unused-batch", {});
        return 9;
    } catch (const std::exception&) {}
    try {
        equilibrium::sample_gfile_boozer("missing.g", 9, 24);
        return 10;
    } catch (const std::exception&) {}
    std::istringstream malformed("invalid gfile");
    try {
        equilibrium::gfile_to_cumes_input(malformed);
        return 8;
    } catch (const std::exception&) {}
#endif
    equilibrium::Spdata spdata;
    if (equilibrium::validate_spdata(spdata)) return 5;
    equilibrium::Equilibrium owner;
    const auto view = equilibrium::as_view(owner);
    if (!view.families.empty()) return 1;
#ifdef EQUILIBRIUM_COMPONENT_binary
    try {
        equilibrium::read_boozer_binary("does-not-exist.bin");
        return 4;
    } catch (const std::exception&) {}
    if (equilibrium::read_binary("does-not-exist.bin")) return 2;
    if (equilibrium::read_spectral_payload("does-not-exist.bin")) return 3;
#endif
#ifdef EQUILIBRIUM_COMPONENT_io
    if (equilibrium::read("does-not-exist.bin")) return 2;
    if (equilibrium::write_spdata("invalid-spdata.dat", spdata)) return 6;
#endif
}
