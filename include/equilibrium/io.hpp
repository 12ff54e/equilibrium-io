#pragma once

#include "equilibrium/binary.hpp"

namespace equilibrium {

// Inspect the container signature; unknown containers never fall back to bin.
// NetCDF-4 and HDF5 share a signature: use .nc or the explicit-format read.
Result<OutputFormat> detect_format(const std::string& path);
Result<Equilibrium> read(const std::string& path,
                         const ReadRequest& request = ReadRequest::all());
Result<Equilibrium> read(const std::string& path,
                         OutputFormat format,
                         const ReadRequest& request = ReadRequest::all());
Status write(const std::string& path,
             const EquilibriumView& view,
             WriteOptions options = {});
inline Status write(const std::string& path,
                    const Equilibrium& owner,
                    WriteOptions options = {}) {
    return write(path, as_view(owner), options);
}
Status convert(const std::string& input, const std::string& output);

}  // namespace equilibrium
