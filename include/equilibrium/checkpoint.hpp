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

#include <string>

namespace equilibrium {

// Write the latest checkpoint schema. Read checkpoints with read_binary().
Status write_checkpoint(const SnapshotView& snapshot,
                        const InputParams& input_params,
                        const std::string& path);

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INCLUDE_IO_CHECKPOINT_HPP_
