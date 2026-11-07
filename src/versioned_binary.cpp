// versioned_binary.cpp — the schema v1 binary state container (blueprint
// §6.13) and the host-library binary factories.
//
// Layout (little-endian), symmetric version 8 (asymmetric extension: version
// 10):
//   magic     8 bytes  "CUMES001"
//   version   int32    = 8
//   ns        int32
//   mnmax     int32
//   families  6 * (mnmax*ns) doubles, mode-major (the state payload)
//   fields    ntheta(i32), nzeta(i32), then 7 half-grid and 6 full-grid
//             real-space double arrays (version 8; 0,0 means absent)
//   ---- provenance trailer (after the state, so a reader can stop early) ----
//   precision int32    (0=double, 1=float; records the computation type)
//   status    int32    (RunStatus)
//   total_iter int32
//   nstages   int32
//   build     revision(str), dirty(u8), build_type(str),
//             precision_policy(str), compile_flags(str)
//   input     source_path(str), source_hash(str)
//   runtime   gpu_name(str), driver(str), runtime(str), toolkit(str)
//   stages    per stage: ns(i32), iterations(i32), converged(u8),
//             fsqr(f64), fsqz(f64), fsql(f64), nrestarts(i32), restarts(i32...)
//   params    the embedded normalized-input record (io_common.hpp
//             write_input_params), the LAST trailer element; version 3 and up
//
// All strings are int32-length-prefixed. The scalar_type string is NOT
// serialized: the reader reconstructs it from the precision tag.
//
// Version 1 (historical, still readable): the trailer instead carried a
// scalar_type string between build_type and source_path (the v1 reader
// consumes and discards it). Versions 1 and 2 carry no input record; the
// reader reports a default-empty InputParams for them. Version 3 carries the
// input record WITHOUT the three profile-type strings; version 4 appends
// pmass_type/piota_type/pcurr_type after the schema tag (the reader keeps
// the "power_series" defaults for version-3 records). Free-boundary versions
// 5 and 6 add inline-Makegrid provenance; version 7 combines both extensions;
// version 8 inserts scientific fields; version 9 adds a family count after
// ns/mnmax. Version 10 stores asymmetric input as typed vectors; the version-9
// JSON input record is no longer supported (its state remains readable).
// Version 11 appends pflux_type and carries a family count for either symmetry.
//
// The state payload is read and validated independently of the provenance
// trailer, so a reader stays forward-compatible with later v1.x trailers.
#include "equilibrium/reader.hpp"
#include "equilibrium/writer.hpp"
#include "io_common.hpp"
#include "read_selection.hpp"

#include <cstdio>
#include <cstring>
#include <memory>
#include <optional>
#include <string>

namespace equilibrium {
namespace {

constexpr char MAGIC[9] = "CUMES001";
// v2 adds the precision-policy provenance fields to the trailer; v1 files
// are still read (the fields default empty). v3 appends the embedded
// normalized-input record as the last trailer element; v1/v2 files are still
// read (the record defaults empty). v4 appends the three profile-type
// strings to that record; v3 files are still read (the types default to
// "power_series"). Versions 5/6 are the free-boundary lineage; version 7
// combines the profile and free-boundary extensions. Version 8 inserts the
// scientific fields after the stable spectral state payload.
constexpr std::int32_t VERSION = 11;

// The on-disk precision discriminator of the v1 trailer (0=double, 1=float).
// Typed locally instead of a string compare against
// BuildProvenance::scalar_type so an unknown tag is a write error rather than a
// silent "0=double" record.
enum class PrecisionTag : std::int32_t { DOUBLE = 0, FLOAT = 1 };

std::optional<PrecisionTag> precision_tag_for(const std::string& scalar_type) {
    if (scalar_type == "double") return PrecisionTag::DOUBLE;
    if (scalar_type == "float") return PrecisionTag::FLOAT;
    return std::nullopt;
}

class VersionedBinaryWriter final : public Writer {
   public:
    Status write_atomic(const SnapshotView& snapshot,
                        const RunReport& report,
                        const OutputSpec& spec) override {
        if (snapshot.lasym() != report.input_params.lasym) {
            return Status("versioned binary: state/input symmetry mismatch");
        }
        const std::string tmp = io_detail::temp_path_for(spec.path);
        FILE* fp = fopen(tmp.c_str(), "wb");
        if (!fp) return Status("cannot open " + tmp + " for writing");

        auto fail = [&](const std::string& reason) -> Status {
            fclose(fp);
            remove(tmp.c_str());
            return Status(reason);
        };

        const bool flux_extension =
            report.input_params.pflux_type != "power_series";
        const int version =
            flux_extension ? VERSION : (snapshot.lasym() ? 10 : 8);
        bool ok =
            io_detail::write_bytes(fp, MAGIC, 8) &&
            io_detail::write_i32(fp, version) &&
            io_detail::write_i32(fp, snapshot.ns) &&
            io_detail::write_i32(fp, snapshot.mnmax) &&
            (version < 9 || io_detail::write_i32(fp, snapshot.components()));
        if (!ok) return fail("failed to write versioned state payload");
        // write_state_families aborts on a family-size mismatch before writing
        // (an undersized family must not fall through into an OOB read).
        if (!io_detail::write_state_families(fp, snapshot)) {
            return fail("failed to write versioned state payload");
        }
        if (!io_detail::write_derived_fields(fp, snapshot)) {
            return fail("failed to write derived field payload");
        }

        const auto precision = precision_tag_for(report.build.scalar_type);
        if (!precision)
            return fail("unknown precision tag '" + report.build.scalar_type +
                        "'");
        ok = io_detail::write_i32(fp, static_cast<std::int32_t>(*precision)) &&
             io_detail::write_i32(fp,
                                  static_cast<std::int32_t>(report.status)) &&
             io_detail::write_i32(fp, report.total_effective_iterations) &&
             io_detail::write_i32(
                 fp, static_cast<std::int32_t>(report.stages.size()));
        ok = ok && io_detail::write_string(fp, report.build.revision) &&
             io_detail::write_u8(fp, report.build.dirty ? 1 : 0) &&
             io_detail::write_string(fp, report.build.build_type) &&
             io_detail::write_string(fp, report.build.precision_policy) &&
             io_detail::write_string(fp, report.build.compile_flags) &&
             io_detail::write_string(fp, report.input.source_path) &&
             io_detail::write_string(fp, report.input.source_hash) &&
             io_detail::write_string(fp, report.runtime.gpu_name) &&
             io_detail::write_string(fp, report.runtime.driver) &&
             io_detail::write_string(fp, report.runtime.runtime) &&
             io_detail::write_string(fp, report.runtime.toolkit);
        for (const auto& stage : report.stages) {
            ok = ok && io_detail::write_i32(fp, stage.ns) &&
                 io_detail::write_i32(fp, stage.effective_iterations) &&
                 io_detail::write_u8(fp, stage.converged ? 1 : 0) &&
                 io_detail::write_f64(fp, stage.final_residual.fsqr) &&
                 io_detail::write_f64(fp, stage.final_residual.fsqz) &&
                 io_detail::write_f64(fp, stage.final_residual.fsql) &&
                 io_detail::write_i32(
                     fp, static_cast<std::int32_t>(stage.restarts.size()));
            for (const auto& r : stage.restarts) {
                ok = ok && io_detail::write_i32(fp, r.iteration);
            }
        }
        // The embedded normalized-input record is the LAST trailer element.
        ok = ok &&
             io_detail::write_input_params(fp, report.input_params,
                                           snapshot.lasym(), flux_extension);
        if (!ok) return fail("failed to write versioned provenance trailer");

        const std::string err = io_detail::publish_atomic(fp, tmp, spec.path);
        if (!err.empty()) return Status("versioned binary publish: " + err);
        return Status();
    }
};

// Read framing even when provenance is ignored: saved grid/input data follow
// it.
bool read_provenance(FILE* fp,
                     int version,
                     RunReport& report,
                     bool keep,
                     std::string& reason) {
    using namespace io_detail;
    std::int32_t precision = 0, status = 0, total = 0, stages = 0;
    if (!read_i32(fp, precision) || !read_i32(fp, status) ||
        !read_i32(fp, total) || !read_i32(fp, stages) || stages < 0 ||
        stages > MAX_INPUT_PARAMS_VECTOR)
        return false;
    if (keep && ((precision != 0 && precision != 1) || status < 0 ||
                 status > 4 || total < 0))
        reason = "invalid provenance tags";
    report.status = static_cast<RunStatus>(status);
    report.total_effective_iterations = total;
    report.build.scalar_type = precision == 0 ? "double" : "float";
    auto text = [&](std::string& value) {
        return keep ? read_string(fp, value) : skip_string(fp);
    };
    std::uint8_t dirty = 0;
    std::string legacy;
    if (!text(report.build.revision) || !read_u8(fp, dirty) ||
        !text(report.build.build_type) ||
        (version >= 2 ? (!text(report.build.precision_policy) ||
                         !text(report.build.compile_flags))
                      : !text(legacy)) ||
        !text(report.input.source_path) || !text(report.input.source_hash) ||
        !text(report.runtime.gpu_name) || !text(report.runtime.driver) ||
        !text(report.runtime.runtime) || !text(report.runtime.toolkit))
        return false;
    if (keep && dirty > 1) reason = "invalid provenance dirty flag";
    report.build.dirty = dirty != 0;
    for (int j = 0; j < stages; ++j) {
        StageReport stage;
        std::uint8_t converged = 0;
        std::int32_t restarts = 0;
        if (!read_i32(fp, stage.ns) ||
            !read_i32(fp, stage.effective_iterations) ||
            !read_u8(fp, converged) ||
            !read_f64(fp, stage.final_residual.fsqr) ||
            !read_f64(fp, stage.final_residual.fsqz) ||
            !read_f64(fp, stage.final_residual.fsql) ||
            !read_i32(fp, restarts) || restarts < 0 ||
            restarts > MAX_INPUT_PARAMS_VECTOR)
            return false;
        if (keep &&
            (stage.ns < 1 || stage.effective_iterations < 0 || converged > 1))
            reason = "invalid provenance stage tags";
        stage.converged = converged != 0;
        for (int k = 0; k < restarts; ++k) {
            std::int32_t iteration = 0;
            if (!read_i32(fp, iteration)) return false;
            if (keep) {
                if (iteration < 0) reason = "negative restart iteration";
                stage.restarts.push_back({iteration});
            }
        }
        if (keep) report.stages.push_back(std::move(stage));
    }
    return true;
}

class VersionedBinaryReader final : public Reader {
   public:
    Result<EquilibriumSnapshot> read(
        const std::string& path,
        std::optional<std::reference_wrapper<RunReport>> report) override {
        auto decoded = io_detail::read_native_binary(
            path, io_detail::legacy_request(report.has_value()), true);
        if (!decoded) return Result<EquilibriumSnapshot>(decoded.error());
        auto& owner = decoded.value();
        if (owner.source.checkpoint)
            return Result<EquilibriumSnapshot>(
                "expected a result, got a checkpoint");
        if (report) {
            RunReport parsed;
            if (owner.input) parsed.input_params = std::move(*owner.input);
            if (owner.provenance)
                static_cast<RunProvenance&>(parsed) =
                    std::move(*owner.provenance);
            report->get() = std::move(parsed);
        }
        return io_detail::take_snapshot(std::move(owner));
    }
    Result<Equilibrium> read(const std::string& path,
                             const ReadRequest& request) override {
        return io_detail::read_native_binary(path, request);
    }
};
}  // namespace

namespace io_detail {

Result<Equilibrium> read_native_binary(const std::string& path,
                                       const ReadRequest& request,
                                       bool strict) {
    const auto close = [](FILE* fp) { fclose(fp); };
    std::unique_ptr<FILE, decltype(close)> file(fopen(path.c_str(), "rb"),
                                                close);
    if (!file) return Result<Equilibrium>("cannot open " + path);
    FILE* fp = file.get();
    try {
        ReadSelection selection(request, strict);
        Equilibrium owner;
        const auto finish = [&]() -> Result<Equilibrium> {
            if (!strict) selection.reconcile(owner);
            selection.finish();
            owner.read_diagnostics = std::move(selection.diagnostics);
            return std::move(owner);
        };
        const auto check = [&](auto group, bool present) {
            if (!selection.wants(group)) return false;
            // Legacy Reader/checkpoint APIs return raw snapshots. Their owner
            // adapters retain the stricter complete-record validation.
            if (strict) return selection.accept(group, present);
            auto status = [&] {
                if constexpr (std::is_same_v<decltype(group), Quantity>)
                    return validate_quantity(as_view(owner), group);
                else
                    return validate_group(as_view(owner), group);
            }();
            if (!present) return selection.accept(group, false);
            return selection.accept(group, present,
                                    status ? "" : status.error());
        };
        char magic[8];
        std::int32_t version = 0, precision = 0, ns = 0, modes = 0,
                     components = 6;
        if (!read_bytes(fp, magic, 8) || !read_i32(fp, version))
            throw std::runtime_error("truncated binary header");
        const bool checkpoint = std::memcmp(magic, "CUMECKP1", 8) == 0;
        if (!checkpoint && std::memcmp(magic, MAGIC, 8) != 0)
            throw std::runtime_error("unsupported native binary signature");
        if (version < 1 || version > (checkpoint ? 9 : VERSION))
            throw std::runtime_error("unsupported native binary version");
        if ((checkpoint && (!read_i32(fp, precision) || precision != 0)) ||
            !read_i32(fp, ns) || !read_i32(fp, modes) ||
            (version >= (checkpoint ? 7 : 9) && !read_i32(fp, components)))
            throw std::runtime_error("invalid native binary header");
        std::size_t count = 0;
        std::string reason;
        if (!check_state_dimensions(fp, ns, modes, count, reason, components))
            throw std::runtime_error(reason);
        owner.ns = ns;
        owner.mnmax = modes;
        owner.symmetry = components == 12;
        owner.source = {OutputFormat::BINARY, version, checkpoint};
        if (selection.wants(DataGroup::GEOMETRY)) {
            owner.families.resize(components);
            bool ok = true;
            for (std::size_t i = 0; i < owner.families.size(); ++i) {
                auto& family = owner.families[i];
                if (!selection.wants(native_family_quantity(i))) {
                    if (!skip_bytes(fp, count * sizeof(double))) {
                        ok = false;
                        break;
                    }
                    continue;
                }
                family.resize(count);
                if (!read_f64_array(fp, family)) {
                    ok = false;
                    break;
                }
            }
            if (!ok) throw std::runtime_error("truncated spectral geometry");
            for (auto q : group_quantities(DataGroup::GEOMETRY))
                if (!check(q, true)) drop_quantity(owner, q);
        } else if (!skip_bytes(fp, count * components * sizeof(double))) {
            throw std::runtime_error(
                "cannot locate groups after spectral geometry");
        }
        const bool metadata = selection.wants(DataGroup::GRID) ||
                              selection.wants(DataGroup::INPUT);
        if (!checkpoint && (metadata || selection.wants(DataGroup::FIELDS) ||
                            selection.wants(DataGroup::PROVENANCE))) {
            EquilibriumSnapshot fields;
            fields.ns = ns;
            reason.clear();
            if (version >= 8 &&
                !read_derived_fields(fp, fields, reason,
                                     selection.wants(DataGroup::FIELDS),
                                     &request)) {
                selection.unavailable(DataGroup::FIELDS, reason);
                for (auto group :
                     {DataGroup::GRID, DataGroup::INPUT, DataGroup::PROVENANCE})
                    selection.unavailable(
                        group, "cannot locate group after fields: " + reason);
                return finish();
            }
            if (fields.has_any_derived_fields())
                owner.fields = {fields.ntheta, fields.nzeta,
                                std::move(fields.half_fields),
                                std::move(fields.full_fields)};
            for (auto q : group_quantities(DataGroup::FIELDS))
                if (!check(q, owner.fields.has_value()))
                    drop_quantity(owner, q);
            if (metadata || selection.wants(DataGroup::PROVENANCE)) {
                RunReport report;
                reason.clear();
                const bool framed = read_provenance(
                    fp, version, report, selection.wants(DataGroup::PROVENANCE),
                    reason);
                if (!framed) {
                    reason = "truncated or corrupt provenance framing";
                    selection.unavailable(DataGroup::PROVENANCE, reason);
                    selection.unavailable(DataGroup::GRID, reason);
                    selection.unavailable(DataGroup::INPUT, reason);
                    return finish();
                }
                if (selection.accept(DataGroup::PROVENANCE, true, reason))
                    owner.provenance =
                        std::move(static_cast<RunProvenance&>(report));
            }
        }
        const bool has_input = version >= (checkpoint ? 2 : 3);
        if (strict && !checkpoint && metadata && !has_input &&
            (fgetc(fp) != EOF || ferror(fp)))
            throw std::runtime_error("trailing data or read error");
        if (metadata && has_input) {
            if (version == (checkpoint ? 7 : 9)) {
                reason = "JSON-based input record is no longer supported";
                selection.unavailable(DataGroup::GRID, reason);
                selection.unavailable(DataGroup::INPUT, reason);
                return finish();
            }
            const long start = ftell(fp);
            GridMetadata grid;
            if (!read_i32(fp, grid.mpol) || !read_i32(fp, grid.ntor) ||
                !read_i32(fp, grid.nfp) || !read_i32(fp, grid.ntheta) ||
                !read_i32(fp, grid.nzeta)) {
                selection.unavailable(DataGroup::GRID,
                                      "truncated grid metadata");
                selection.unavailable(DataGroup::INPUT,
                                      "truncated input metadata");
                return finish();
            }
            owner.grid = grid;
            if (!check(DataGroup::GRID, true)) owner.grid.reset();
            if (selection.wants(DataGroup::INPUT)) {
                InputParams input;
                reason.clear();
                const bool profile_types = checkpoint
                                               ? version == 3 || version >= 6
                                               : version == 4 || version >= 7;
                const int free_boundary = checkpoint
                                              ? (version >= 6   ? 3
                                                 : version >= 4 ? version - 2
                                                                : 0)
                                              : (version >= 7   ? 3
                                                 : version >= 5 ? version - 3
                                                                : 0);
                const bool asymmetric =
                    checkpoint
                        ? version == 8 || (version >= 9 && components == 12)
                        : version == 10 || (version >= 11 && components == 12);
                const bool decoded =
                    start >= 0 && fseek(fp, start, SEEK_SET) == 0 &&
                    read_input_params(fp, input, reason, profile_types,
                                      free_boundary, asymmetric,
                                      version >= (checkpoint ? 9 : 11));
                if (!decoded) {
                    selection.unavailable(
                        DataGroup::INPUT,
                        reason.empty() ? "truncated input record" : reason);
                } else {
                    owner.input = std::move(input);
                    if (strict && owner.input->lasym != (components == 12)) {
                        selection.unavailable(DataGroup::INPUT,
                                              "state/input symmetry mismatch");
                        owner.input.reset();
                    } else if (!check(DataGroup::INPUT, true)) {
                        owner.input.reset();
                    }
                    if (strict && (fgetc(fp) != EOF || ferror(fp)))
                        throw std::runtime_error("trailing data or read error");
                }
            }
        }
        return finish();
    } catch (const std::exception& error) {
        return Result<Equilibrium>("native binary: " +
                                   std::string(error.what()));
    }
}

}  // namespace io_detail

// Binary-only factories are available without optional file libraries.
std::unique_ptr<Writer> make_binary_writer() {
    return std::make_unique<VersionedBinaryWriter>();
}

std::unique_ptr<Reader> make_binary_reader() {
    return std::make_unique<VersionedBinaryReader>();
}

}  // namespace equilibrium
