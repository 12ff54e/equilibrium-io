#include "equilibrium/checkpoint.hpp"
#include "equilibrium/io.hpp"
#include "equilibrium/validation.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>

namespace eq = equilibrium;

template <class T>
concept BorrowableTemporary = requires { eq::as_view(T{}); };
static_assert(!BorrowableTemporary<eq::Equilibrium>);
static_assert(!BorrowableTemporary<eq::EquilibriumSnapshot>);
static_assert(!BorrowableTemporary<eq::EquilibriumProfiles>);

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

eq::Equilibrium fixture() {
    eq::Equilibrium owner;
    owner.ns = 3;
    owner.mnmax = 2;
    owner.families.resize(6);
    for (std::size_t c = 0; c < owner.families.size(); ++c) {
        owner.families[c].resize(6);
        for (std::size_t i = 0; i < 6; ++i)
            owner.families[c][i] = 0.25 * (1 + c * 6 + i);
    }
    auto& input = owner.input.emplace();
    input.mpol = 2;
    input.nfp = 3;
    input.ntheta = 8;
    input.nzeta = 1;
    input.delt = 0.9;
    input.phiedge = input.bloat = input.pres_scale = 1;
    input.nvacskip = 1;
    input.am.assign(1, 0);
    input.ai.assign(1, 0.4);
    input.pcurr_type = "two_power";
    input.ac = {1, 2, 3};
    input.aphi.assign(1, 1);
    input.raxis_c.assign(1, 10);
    input.zaxis_s.assign(1, 0);
    input.rbc_m = {0, 1, 1};
    input.rbc_n = {0, 0, 0};
    input.rbc_value = {10, 0.5, 0};  // Preserve explicit zeros and duplicates.
    input.zbs_m.assign(1, 1);
    input.zbs_n.assign(1, 0);
    input.zbs_value.assign(1, 0.5);
    input.rbcc = {10, 0.5};
    input.rbss = {0, 0};
    input.zbsc = {0, 0.5};
    input.zbcs = {0, 0};
    input.stages.push_back({3, 20, 1e-8});
    owner.grid = eq::grid_metadata(input);
    auto& provenance = owner.provenance.emplace();
    provenance.build.scalar_type = "double";
    provenance.build.revision = "standalone-fixture";
    provenance.status = eq::RunStatus::CONVERGED;
    provenance.stages.push_back({3, 7, true, {1e-9, 2e-9, 3e-9}, {{2}}});
    provenance.total_effective_iterations = 7;
    return owner;
}

std::string bytes(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(file), {}};
}

int main(int argc, char** argv) {
    try {
        require(argc == 2, "expected scratch directory");
        const auto scratch = std::filesystem::path(argv[1]);
        std::filesystem::create_directories(scratch);
        auto owner = fixture();
        require(bool(eq::validate(eq::as_view(owner))),
                "valid fixture rejected");
        auto view = eq::as_view(owner);
        require(view.families[0].data() == owner.families[0].data(),
                "view copied state");
        require(!view.fields && !view.profiles,
                "missing groups were fabricated");
        require(&view.input->get() == &*owner.input,
                "view copied input metadata");
        auto copy = owner;
        require(eq::as_view(copy).families[0].data() != view.families[0].data(),
                "copied owner retained source spans");
        auto moved = std::move(copy);
        require(eq::as_view(moved).families[0][2] == owner.families[0][2],
                "moved owner lost state");

        auto& fields = owner.fields.emplace();
        fields.ntheta = 8;
        fields.nzeta = 1;
        for (auto& field : fields.half) field.assign(16, 0);
        for (auto& field : fields.full) field.assign(24, 0);
        require(bool(eq::validate(eq::as_view(owner))),
                "zero fields are valid data");
        fields.half[0].pop_back();
        require(!eq::validate(eq::as_view(owner)), "partial fields accepted");
        fields.half[0].push_back(0);
        auto& profiles = owner.profiles.emplace();
        profiles.toroidal_flux_derivative = {1, 1};
        profiles.poloidal_flux_derivative = {.4, .4};
        profiles.rotational_transform = {.4, .4};
        profiles.poloidal_covariant_field = {0, 0};
        profiles.toroidal_covariant_field = {1, 1};
        require(eq::as_view(owner).profiles->rotational_transform.data() ==
                    profiles.rotational_transform.data(),
                "view copied profiles");
        const auto binary = (scratch / "saved.bin").string();
        require(!eq::write_binary(binary, owner),
                "profiles silently discarded");
        require(bool(eq::write_binary(binary, owner,
                                      {eq::ProfileWritePolicy::OMIT})),
                "explicit profile omission failed");
        owner.profiles.reset();
        auto saved = eq::read(binary);
        require(bool(saved), "complete binary read failed");
        require(!saved.value().profiles && saved.value().fields.has_value(),
                "decoded field presence changed");
        require(saved.value().source.version == 8,
                "wrong binary source version");
        auto state = eq::read(binary, eq::ReadMode::STATE_ONLY);
        require(bool(state) && !state.value().grid && !state.value().input &&
                    !state.value().provenance,
                "state-only read invented metadata");
        const auto original = bytes(binary);
        owner.input->rbcc[0] += 1;
        require(!eq::write(binary, owner),
                "inconsistent folded boundary accepted");
        require(bytes(binary) == original, "failed write replaced destination");
        owner.input->rbcc[0] -= 1;

        for (auto format : {eq::OutputFormat::BINARY, eq::OutputFormat::NETCDF,
                            eq::OutputFormat::HDF5}) {
            if (!eq::output_format_available(format)) continue;
            const auto path =
                (scratch /
                 ("roundtrip" + std::string(eq::output_suffix(format))))
                    .string();
            const auto status = eq::write(path, owner);
            if (!status) throw std::runtime_error(status.error());
            auto decoded = eq::read(path);
            if (!decoded) throw std::runtime_error(decoded.error());
            require(decoded.value().families == owner.families,
                    "state changed across containers");
            require(decoded.value().fields->half == owner.fields->half,
                    "fields changed");
            require(decoded.value().input == owner.input,
                    "input provenance changed");
            require(decoded.value().provenance->build.revision ==
                        "standalone-fixture",
                    "report lost");
            const auto back = (scratch / "back.bin").string();
            require(bool(eq::convert(path, back)), "conversion failed");
            require(bytes(back) == original, "binary conversion changed bytes");
        }
        const auto checkpoint = (scratch / "checkpoint.bin").string();
        require(
            bool(eq::write_checkpoint(eq::as_snapshot_view(eq::as_view(owner)),
                                      *owner.input, checkpoint)),
            "checkpoint write failed");
        auto restarted = eq::read_binary(checkpoint);
        require(bool(restarted) && restarted.value().source.checkpoint &&
                    !restarted.value().provenance &&
                    !restarted.value().fields &&
                    restarted.value().input == owner.input,
                "checkpoint presence contract broken");
        require(!eq::write(binary, restarted.value()),
                "checkpoint fabricated run provenance");
        const auto unknown = (scratch / "unknown.bin").string();
        std::ofstream(unknown) << "arbitrary data";
        require(!eq::read(unknown), "unknown magic accepted");
        std::ofstream(unknown, std::ios::binary)
            << original.substr(0, original.size() - 3);
        require(!eq::read(unknown), "truncated provenance accepted");
        require(bool(eq::read(unknown, eq::ReadMode::STATE_ONLY)),
                "explicit state-only legacy read failed");
        std::ofstream(unknown, std::ios::binary) << original << 'x';
        require(!eq::read(unknown), "trailing garbage accepted");
        // Header, six spectral families, field dimensions, seven half-grid
        // fields and six full-grid fields precede the provenance tags.
        const std::size_t trailer = 20 + 6 * 6 * 8 + 8 + (7 * 16 + 6 * 24) * 8;
        for (auto offset : {trailer, trailer + 4}) {
            auto corrupt = original;
            corrupt[offset] = 0;
            corrupt[offset + 1] = 1;  // 256 must not narrow to CONVERGED.
            std::ofstream(unknown, std::ios::binary)
                .write(corrupt.data(), corrupt.size());
            require(!eq::read(unknown), "invalid binary tag accepted");
            require(bool(eq::read(unknown, eq::ReadMode::STATE_ONLY)),
                    "invalid provenance prevented state-only read");
        }
        std::cout << "equilibrium contracts passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
