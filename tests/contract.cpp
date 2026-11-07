#include "equilibrium/checkpoint.hpp"
#include "equilibrium/io.hpp"
#include "equilibrium/validation.hpp"
#include "equilibrium/writer.hpp"

#include <cstdint>
#ifdef EQUILIBRIUM_HAVE_NETCDF
#include <netcdf.h>
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
#include <hdf5.h>
#endif

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

void test_requests(const std::filesystem::path& scratch) {
    using G = eq::DataGroup;
    using R = eq::Requirement;
    using A = eq::Availability;
    auto owner = fixture();
    auto& fields = owner.fields.emplace();
    fields.ntheta = 8;
    fields.nzeta = 1;
    for (auto& field : fields.half) field.assign(16, 1);
    for (auto& field : fields.full) field.assign(24, 2);
    eq::ReadRequest request;
    for (auto group : {G::GEOMETRY, G::FIELDS, G::GRID})
        request[group] = R::REQUIRED;
    for (auto format : {eq::OutputFormat::BINARY, eq::OutputFormat::NETCDF,
                        eq::OutputFormat::HDF5}) {
        if (!eq::output_format_available(format)) continue;
        const auto path =
            (scratch / ("selected" + std::string(eq::output_suffix(format))))
                .string();
        const auto reset = [&] {
            auto status = eq::write(path, owner);
            if (!status) throw std::runtime_error(status.error());
        };
        reset();
        auto decoded = eq::read(path, request);
        require(bool(decoded), "requested physical groups failed");
        require(
            decoded.value().has_geometry() && decoded.value().has_fields() &&
                decoded.value().has_grid() && !decoded.value().has_input() &&
                !decoded.value().has_provenance(),
            "ignored groups exposed");
        require(decoded.value().read_diagnostics[G::INPUT].availability ==
                    A::NOT_REQUESTED,
                "ignored status lost");
        const auto corrupt = [&](bool provenance) {
            if (format == eq::OutputFormat::BINARY) {
                auto data = bytes(path);
                const auto offset =
                    provenance ? 20 + 6 * 6 * 8 + 8 + (7 * 16 + 6 * 24) * 8
                               : data.find("cumes-config-v1") - 72;
                data[offset] = 7;
                data[offset + 1] = data[offset + 2] = data[offset + 3] = 0;
                std::ofstream(path, std::ios::binary)
                    .write(data.data(), data.size());
            }
#ifdef EQUILIBRIUM_HAVE_NETCDF
            if (format == eq::OutputFormat::NETCDF) {
                int file = -1, var = -1, value = 7;
                require(nc_open(path.c_str(), NC_WRITE, &file) == NC_NOERR,
                        "nc open");
                require(nc_inq_varid(file, provenance ? "precision" : "ncurr",
                                     &var) == NC_NOERR,
                        "nc var");
                require(nc_put_var_int(file, var, &value) == NC_NOERR,
                        "nc corrupt");
                require(nc_close(file) == NC_NOERR, "nc close");
            }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
            if (format == eq::OutputFormat::HDF5) {
                auto file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                auto attr = H5Aopen(file, provenance ? "precision" : "ncurr",
                                    H5P_DEFAULT);
                int value = 7;
                require(file >= 0 && attr >= 0 &&
                            H5Awrite(attr, H5T_NATIVE_INT, &value) >= 0,
                        "h5 corrupt");
                H5Aclose(attr);
                H5Fclose(file);
            }
#endif
        };
        for (bool provenance : {true, false}) {
            const auto group = provenance ? G::PROVENANCE : G::INPUT;
            reset();
            corrupt(provenance);
            auto selected = request;
            decoded = eq::read(path, selected);
            require(bool(decoded) && decoded.value().has_grid(),
                    "ignored metadata blocked required grid");
            selected[group] = R::OPTIONAL;
            decoded = eq::read(path, selected);
            require(bool(decoded),
                    "bad optional metadata blocked required groups");
            const auto& status = decoded.value().read_diagnostics[group];
            require(
                status.availability == A::INVALID && !status.diagnostic.empty(),
                "bad optional metadata lost its diagnostic");
            require(!decoded.value().has_input() &&
                        !decoded.value().has_provenance(),
                    "partial invalid metadata exposed");
            selected[group] = R::REQUIRED;
            require(!eq::read(path, selected),
                    "bad required metadata accepted");
        }

        if (format != eq::OutputFormat::BINARY) {
            for (bool provenance : {true, false})
                for (bool remove_identifiers : {false, true}) {
                    reset();
                    std::vector<const char*> names{provenance ? "precision"
                                                              : "mpol"};
                    if (remove_identifiers)
                        names =
                            provenance
                                ? std::vector<const char*>{"precision",
                                                           "status",
                                                           "total_iterations",
                                                           "build_dirty",
                                                           "revision",
                                                           "stage_ns"}
                                : std::vector<const char*>{
                                      "mpol", "ntor", "nfp", "ntheta", "nzeta"};
#ifdef EQUILIBRIUM_HAVE_NETCDF
                    if (format == eq::OutputFormat::NETCDF) {
                        int file = -1;
                        require(nc_open(path.c_str(), NC_WRITE, &file) ==
                                        NC_NOERR &&
                                    nc_redef(file) == NC_NOERR,
                                "nc edit metadata");
                        for (const char* name : names) {
                            int var = -1;
                            if (nc_inq_varid(file, name, &var) == NC_NOERR)
                                require(nc_rename_var(
                                            file, var,
                                            (std::string("removed_") + name)
                                                .c_str()) == NC_NOERR,
                                        "nc remove metadata variable");
                            else
                                require(nc_del_att(file, NC_GLOBAL, name) ==
                                            NC_NOERR,
                                        "nc remove metadata attribute");
                        }
                        require(nc_close(file) == NC_NOERR, "nc close");
                    }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
                    if (format == eq::OutputFormat::HDF5) {
                        auto file =
                            H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                        require(file >= 0, "h5 edit metadata");
                        for (const char* name : names)
                            require(
                                (H5Aexists(file, name) > 0
                                     ? H5Adelete(file, name)
                                     : H5Ldelete(file, name, H5P_DEFAULT)) >= 0,
                                "h5 remove metadata");
                        H5Fclose(file);
                    }
#endif
                    eq::ReadRequest selected;
                    const auto group = provenance ? G::PROVENANCE : G::INPUT;
                    selected[group] = R::OPTIONAL;
                    decoded = eq::read(path, selected);
                    require(
                        bool(decoded) && decoded.value()
                                                 .read_diagnostics[group]
                                                 .availability == A::INVALID,
                        "partially present metadata reported absent");
                }
        }
        if (format != eq::OutputFormat::BINARY) {
            reset();
#ifdef EQUILIBRIUM_HAVE_NETCDF
            if (format == eq::OutputFormat::NETCDF) {
                int file = -1, var = -1;
                require(nc_open(path.c_str(), NC_WRITE, &file) == NC_NOERR,
                        "nc open");
                require(nc_inq_varid(file, "sqrtg", &var) == NC_NOERR,
                        "nc field");
                require(
                    nc_redef(file) == NC_NOERR &&
                        nc_rename_var(file, var, "removed_sqrtg") == NC_NOERR,
                    "nc remove field");
                require(nc_close(file) == NC_NOERR, "nc close");
            }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
            if (format == eq::OutputFormat::HDF5) {
                auto file = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                require(file >= 0 && H5Ldelete(file, "sqrtg", H5P_DEFAULT) >= 0,
                        "h5 remove field");
                H5Fclose(file);
            }
#endif
            auto selected = request;
            selected[G::FIELDS] = R::OPTIONAL;
            decoded = eq::read(path, selected);
            require(
                bool(decoded) && !decoded.value().has_fields() &&
                    decoded.value().has_grid() &&
                    decoded.value().read_diagnostics[G::FIELDS].availability ==
                        A::INVALID,
                "partial optional scientific datasets blocked grid");
            selected[G::FIELDS] = R::IGNORED;
            require(bool(eq::read(path, selected)),
                    "ignored scientific datasets decoded");
            selected[G::FIELDS] = R::REQUIRED;
            require(!eq::read(path, selected),
                    "partial required datasets accepted");
        }

        if (format != eq::OutputFormat::BINARY) {
            auto asymmetric = owner;
            asymmetric.families.resize(12, std::vector<double>(6, 0));
            auto& input = *asymmetric.input;
            input.lasym = true;
            input.raxis_s = input.zaxis_c = {0};
            input.rbsc = input.rbcs = input.zbcc = input.zbss = {0, 0};
            require(bool(eq::write(path, asymmetric)),
                    "asymmetric selection fixture");
            const auto corrupt_asymmetry = [&](bool tag, int value = 2) {
#ifdef EQUILIBRIUM_HAVE_NETCDF
                if (format == eq::OutputFormat::NETCDF) {
                    int file = -1, var = -1;
                    require(nc_open(path.c_str(), NC_WRITE, &file) == NC_NOERR,
                            "nc open");
                    require(nc_inq_varid(file, tag ? "lasym" : "rmnsc", &var) ==
                                NC_NOERR,
                            "nc var");
                    if (tag)
                        require(nc_put_var_int(file, var, &value) == NC_NOERR,
                                "nc tag");
                    else
                        require(nc_redef(file) == NC_NOERR &&
                                    nc_rename_var(file, var, "removed_rmnsc") ==
                                        NC_NOERR,
                                "nc malformed geometry");
                    require(nc_close(file) == NC_NOERR, "nc close");
                }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
                if (format == eq::OutputFormat::HDF5) {
                    auto file =
                        H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    require(file >= 0, "h5 open");
                    if (tag) {
                        auto attr = H5Aopen(file, "lasym", H5P_DEFAULT);
                        require(attr >= 0 &&
                                    H5Awrite(attr, H5T_NATIVE_INT, &value) >= 0,
                                "h5 tag");
                        H5Aclose(attr);
                    } else
                        require(H5Ldelete(file, "rmnsc", H5P_DEFAULT) >= 0,
                                "h5 malformed geometry");
                    H5Fclose(file);
                }
#endif
            };
            corrupt_asymmetry(true, 0);
            eq::ReadRequest symmetric_input;
            symmetric_input[G::GEOMETRY] = R::OPTIONAL;
            symmetric_input[G::INPUT] = R::REQUIRED;
            decoded = eq::read(path, symmetric_input);
            if (!decoded) throw std::runtime_error(decoded.error());
            require(bool(decoded) && !decoded.value().has_geometry() &&
                        decoded.value().has_input() &&
                        !decoded.value().input->lasym,
                    "optional asymmetric geometry invalidated symmetric input");
            require(bool(eq::write(path, asymmetric)),
                    "restore asymmetric fixture");
            corrupt_asymmetry(false);
            eq::ReadRequest selected;
            selected[G::GEOMETRY] = R::OPTIONAL;
            selected[G::INPUT] = R::REQUIRED;
            decoded = eq::read(path, selected);
            require(bool(decoded) && !decoded.value().has_geometry() &&
                        decoded.value().has_input() &&
                        decoded.value().input->lasym &&
                        decoded.value().input->rbsc == input.rbsc,
                    "optional geometry failure lost asymmetric saved input");
            corrupt_asymmetry(true);
            selected[G::GEOMETRY] = R::IGNORED;
            require(!eq::read(path, selected),
                    "malformed input symmetry accepted");
            selected[G::INPUT] = R::OPTIONAL;
            decoded = eq::read(path, selected);
            require(
                bool(decoded) && !decoded.value().has_input() &&
                    decoded.value().read_diagnostics[G::INPUT].availability ==
                        A::INVALID,
                "invalid optional input symmetry lost diagnostic");
        }
        reset();
        auto input_only = eq::ReadRequest{};
        input_only[G::INPUT] = R::REQUIRED;
        decoded = eq::read(path, input_only);
        require(bool(decoded) && decoded.value().has_input() &&
                    !decoded.value().has_geometry() &&
                    !decoded.value().has_fields() &&
                    !decoded.value().has_grid(),
                "input requires exposed unrelated groups");

        // Both groups are locally valid, but the scientific grid and saved
        // angular metadata disagree. Required data take precedence.
        eq::RunReport report;
        static_cast<eq::RunProvenance&>(report) = *owner.provenance;
        report.input_params = *owner.input;
        report.input_params.ntheta = 12;
        require(bool(eq::create_writer(format)->write_atomic(
                    eq::as_snapshot_view(eq::as_view(owner)), report,
                    {format, path})),
                "conflicting group fixture");
        auto conflict = request;
        conflict[G::FIELDS] = R::OPTIONAL;
        decoded = eq::read(path, conflict);
        require(bool(decoded) && decoded.value().has_grid() &&
                    decoded.value().grid->ntheta == 12 &&
                    !decoded.value().has_fields() &&
                    decoded.value().read_diagnostics[G::FIELDS].availability ==
                        A::INVALID,
                "optional fields invalidated required grid");
        conflict[G::FIELDS] = R::REQUIRED;
        conflict[G::GRID] = R::OPTIONAL;
        decoded = eq::read(path, conflict);
        require(bool(decoded) && decoded.value().has_fields() &&
                    !decoded.value().has_grid(),
                "optional grid invalidated required fields");
        conflict[G::GRID] = R::REQUIRED;
        require(!eq::read(path, conflict),
                "conflicting required groups accepted");
        conflict[G::GRID] = R::OPTIONAL;
        conflict[G::FIELDS] = R::OPTIONAL;
        decoded = eq::read(path, conflict);
        require(bool(decoded) && !decoded.value().has_grid() &&
                    !decoded.value().has_fields(),
                "conflicting optional groups exposed");
        conflict[G::GRID] = R::IGNORED;
        conflict[G::INPUT] = R::REQUIRED;
        decoded = eq::read(path, conflict);
        require(bool(decoded) && decoded.value().has_input() &&
                    !decoded.value().has_fields(),
                "optional fields invalidated required input");
        conflict[G::INPUT] = R::OPTIONAL;
        conflict[G::FIELDS] = R::REQUIRED;
        decoded = eq::read(path, conflict);
        require(bool(decoded) && decoded.value().has_fields() &&
                    !decoded.value().has_input(),
                "optional input invalidated required fields");
        auto without_fields = owner;
        without_fields.fields.reset();
        require(bool(eq::write(path, without_fields)), "write absent fields");
        auto selected = request;
        selected[G::FIELDS] = R::OPTIONAL;
        decoded = eq::read(path, selected);
        require(bool(decoded) && !decoded.value().has_fields() &&
                    decoded.value().read_diagnostics[G::FIELDS].availability ==
                        A::ABSENT,
                "absent fields reported incorrectly");
        selected[G::FIELDS] = R::REQUIRED;
        require(!eq::read(path, selected), "missing required fields accepted");
    }
    const auto path = (scratch / "selected.bin").string();
    require(bool(eq::write(path, owner)), "write binary request fixture");
    auto original = bytes(path);
    auto corrupt = original;
    // First half-grid scalar; encode a quiet NaN independently of host byte
    // order.
    const std::size_t offset = 20 + 6 * 6 * 8 + 8;
    for (int j = 0; j < 8; ++j)
        corrupt[offset + j] =
            static_cast<char>(UINT64_C(0x7ff8000000000000) >> (8 * j));
    std::ofstream(path, std::ios::binary).write(corrupt.data(), corrupt.size());
    auto selected = request;
    selected[G::FIELDS] = R::OPTIONAL;
    auto decoded = eq::read_binary(path, selected);
    require(bool(decoded) && !decoded.value().has_fields() &&
                decoded.value().read_diagnostics[G::FIELDS].availability ==
                    A::INVALID,
            "invalid optional physical fields exposed");
    selected[G::FIELDS] = R::IGNORED;
    decoded = eq::read_binary(path, selected);
    require(bool(decoded) && decoded.value().has_grid(),
            "ignored fields validated");
    selected[G::FIELDS] = R::REQUIRED;
    require(!eq::read_binary(path, selected),
            "nonfinite required fields accepted");
    // A truncated optional field block can be diagnosed when no later group is
    // needed.
    std::ofstream(path, std::ios::binary).write(original.data(), offset + 5);
    selected[G::FIELDS] = R::OPTIONAL;
    selected[G::GRID] = R::IGNORED;
    decoded = eq::read_binary(path, selected);
    require(bool(decoded) &&
                decoded.value().read_diagnostics[G::FIELDS].availability ==
                    A::INVALID,
            "truncated optional field block lost status");
    selected[G::GRID] = R::REQUIRED;
    require(!eq::read_binary(path, selected),
            "unlocatable required grid accepted");

    const auto checkpoint = (scratch / "selected-checkpoint.bin").string();
    require(bool(eq::write_checkpoint(eq::as_snapshot_view(eq::as_view(owner)),
                                      *owner.input, checkpoint)),
            "checkpoint fixture");
    selected = request;
    selected[G::FIELDS] = R::OPTIONAL;
    decoded = eq::read_binary(checkpoint, selected);
    require(bool(decoded) && decoded.value().has_grid() &&
                decoded.value().read_diagnostics[G::FIELDS].availability ==
                    A::ABSENT,
            "checkpoint capability request failed");
}

int main(int argc, char** argv) {
    try {
        require(argc == 2, "expected scratch directory");
        const auto scratch = std::filesystem::path(argv[1]);
        std::filesystem::create_directories(scratch);
        test_requests(scratch);
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
