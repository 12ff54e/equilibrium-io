// checkpoint.cpp — versioned checkpoint.
//
// Versioned checkpoint layout (little-endian):
//   magic     8 bytes  "CUMECKP1"
//   version   int32    = 6
//   precision int32    (0 = double; the checkpoint is always double on disk)
//   ns        int32
//   mnmax     int32
//   families  6 * (mnmax*ns) doubles, mode-major
//   params    the embedded normalized-input record (io_common.hpp
//             write_input_params); version 2 and up
//
// Version 3 appends the three profile-type strings to the input record;
// version-2 checkpoints are still read (the types default to
// "power_series"). Free-boundary versions 4/5 add inline-Makegrid provenance;
// version 6 combines both extensions. The restart path reads the state only and
// ignores the record; version-1 checkpoints (no record) remain readable.
// Version 7 adds the family count; version 8 stores complementary input as
// typed vectors. The version-7 JSON input record is no longer supported.
// Version 9 appends pflux_type and carries a family count for either symmetry.
#include "equilibrium/checkpoint.hpp"

#include "io_common.hpp"
#include "read_selection.hpp"

#include <cstdio>
#include <cstring>
#include <string>

namespace equilibrium {
namespace {

constexpr char CHECKPOINT_MAGIC[9] = "CUMECKP1";
constexpr std::int32_t CHECKPOINT_VERSION = 9;

}  // namespace

Status write_checkpoint(const SnapshotView& snapshot,
                        const InputParams& input_params,
                        const std::string& path) {
    if (snapshot.lasym() != input_params.lasym) {
        return Status("checkpoint: state/input symmetry mismatch");
    }
    const std::string tmp = io_detail::temp_path_for(path);
    FILE* fp = fopen(tmp.c_str(), "wb");
    if (!fp) return Status("cannot open " + tmp + " for writing");

    auto fail = [&](const std::string& reason) -> Status {
        fclose(fp);
        remove(tmp.c_str());
        return Status(reason);
    };

    const bool flux_extension = input_params.pflux_type != "power_series";
    const int version =
        flux_extension ? CHECKPOINT_VERSION : (snapshot.lasym() ? 8 : 6);
    bool ok = io_detail::write_bytes(fp, CHECKPOINT_MAGIC, 8) &&
              io_detail::write_i32(fp, version) &&
              io_detail::write_i32(fp, 0 /* precision = double */) &&
              io_detail::write_i32(fp, snapshot.ns) &&
              io_detail::write_i32(fp, snapshot.mnmax) &&
              (version < 7 || io_detail::write_i32(fp, snapshot.components()));
    if (!ok) return fail("failed to write checkpoint payload");
    // write_state_families aborts on a family-size mismatch before writing (an
    // undersized family must not fall through into an OOB read).
    if (!io_detail::write_state_families(fp, snapshot)) {
        return fail("failed to write checkpoint payload");
    }
    if (!io_detail::write_input_params(fp, input_params, snapshot.lasym(),
                                       flux_extension)) {
        return fail("failed to write checkpoint input record");
    }

    const std::string err = io_detail::publish_atomic(fp, tmp, path);
    if (!err.empty()) return Status("checkpoint publish: " + err);
    return Status();
}

Result<EquilibriumSnapshot> read_checkpoint(
    const std::string& path,
    std::optional<std::reference_wrapper<InputParams>> input_params) {
    auto decoded = io_detail::read_native_binary(
        path, io_detail::legacy_request(input_params.has_value()), true);
    if (!decoded) return Result<EquilibriumSnapshot>(decoded.error());
    auto& owner = decoded.value();
    if (!owner.source.checkpoint)
        return Result<EquilibriumSnapshot>("expected a checkpoint");
    if (input_params)
        input_params->get() =
            owner.input ? std::move(*owner.input) : InputParams{};
    return io_detail::take_snapshot(std::move(owner));
}

}  // namespace equilibrium
