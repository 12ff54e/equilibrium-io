// Checkpoint v9: header, family count, state and typed saved input.
#include "equilibrium/checkpoint.hpp"

#include "format_versions.hpp"
#include "io_common.hpp"

#include <cstdio>
#include <cstring>
#include <string>

namespace equilibrium {
namespace {

constexpr char CHECKPOINT_MAGIC[9] = "CUMECKP1";

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

    bool ok = io_detail::write_bytes(fp, CHECKPOINT_MAGIC, 8) &&
              io_detail::write_i32(fp, io_detail::CHECKPOINT_VERSION) &&
              io_detail::write_i32(fp, 0 /* precision = double */) &&
              io_detail::write_i32(fp, snapshot.ns) &&
              io_detail::write_i32(fp, snapshot.mnmax) &&
              io_detail::write_i32(fp, snapshot.components());
    if (!ok) return fail("failed to write checkpoint payload");
    // write_state_families aborts on a family-size mismatch before writing (an
    // undersized family must not fall through into an OOB read).
    if (!io_detail::write_state_families(fp, snapshot)) {
        return fail("failed to write checkpoint payload");
    }
    if (!io_detail::write_input_params(fp, input_params, snapshot.lasym())) {
        return fail("failed to write checkpoint input record");
    }

    const std::string err = io_detail::publish_atomic(fp, tmp, path);
    if (!err.empty()) return Status("checkpoint publish: " + err);
    return Status();
}

}  // namespace equilibrium
