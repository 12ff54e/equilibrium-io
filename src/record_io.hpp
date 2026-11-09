#pragma once

#include "equilibrium/binary.hpp"
#include "equilibrium/reader.hpp"
#include "equilibrium/writer.hpp"

namespace equilibrium::io_detail {

Status write_record(Writer& writer,
                    const OutputSpec& output,
                    const EquilibriumView& view,
                    WriteOptions options);

}  // namespace equilibrium::io_detail
