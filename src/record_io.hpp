#pragma once

#include "equilibrium/binary.hpp"
#include "equilibrium/reader.hpp"
#include "equilibrium/writer.hpp"

namespace equilibrium::io_detail {

Result<Equilibrium> read_record(Reader& reader,
                                const std::string& path,
                                OutputFormat format,
                                ReadMode mode);
Status write_record(Writer& writer,
                    const OutputSpec& output,
                    const EquilibriumView& view,
                    WriteOptions options);

}  // namespace equilibrium::io_detail
