// writer.hpp — host-memory writer interface (blueprint §6.13).
//
// Writers consume host memory only (no CUDA calls) and publish atomically:
// write to a same-directory temp file, flush + fsync + close, then rename over
// the target. A reader never sees a half-written file. Pre-rename failures
// preserve the target; a post-rename sync error reports unconfirmed durability.
#ifndef EQUILIBRIUM_INCLUDE_IO_WRITER_HPP_
#define EQUILIBRIUM_INCLUDE_IO_WRITER_HPP_

#include "equilibrium/equilibrium_snapshot.hpp"
#include "equilibrium/output_spec.hpp"
#include "equilibrium/result.hpp"
#include "equilibrium/run_report.hpp"

#include <memory>

namespace equilibrium {

class Writer {
   public:
    virtual ~Writer() = default;
    // Writers consume host memory only (no CUDA calls) and publish atomically:
    // write to a same-directory temp file, flush + fsync + close, then rename
    // over the target. The report supplies saved input and provenance.
    virtual Status write_atomic(const SnapshotView& snapshot,
                                const RunReport& report,
                                const OutputSpec& spec) = 0;
};

// Full codec dispatch, linked through equilibrium::io. An unavailable codec
// returns nullptr. Use make_binary_writer with equilibrium::binary alone.
std::unique_ptr<Writer> create_writer(OutputFormat format);

// Binary-only factory (host I/O library; no backend defines required).
std::unique_ptr<Writer> make_binary_writer();

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INCLUDE_IO_WRITER_HPP_
