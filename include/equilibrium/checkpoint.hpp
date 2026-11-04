// checkpoint.hpp — versioned checkpoint.
//
// An explicit, versioned restart artifact: a self-describing checkpoint
// (magic, version, precision, dimensions) that a reader validates before
// accepting.
#ifndef EQUILIBRIUM_INCLUDE_IO_CHECKPOINT_HPP_
#define EQUILIBRIUM_INCLUDE_IO_CHECKPOINT_HPP_

#include "equilibrium/equilibrium_snapshot.hpp"
#include "equilibrium/input_params.hpp"
#include "equilibrium/result.hpp"

#include <functional>
#include <optional>
#include <string>

namespace equilibrium {

// Write a versioned (schema v1) checkpoint. Version 2 appends the embedded
// normalized-input record after the state families.
Status write_checkpoint(const SnapshotView& snapshot,
                        const InputParams& input_params,
                        const std::string& path);

// Read a versioned checkpoint, validating the magic, version, and dimensions.
// The restart path consumes the state only; when `input_params` is engaged
// and the checkpoint carries the version-2 record, it is decoded into it
// (left default-empty for version-1 checkpoints).
Result<EquilibriumSnapshot> read_checkpoint(
    const std::string& path,
    std::optional<std::reference_wrapper<InputParams>> input_params =
        std::nullopt);

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INCLUDE_IO_CHECKPOINT_HPP_
