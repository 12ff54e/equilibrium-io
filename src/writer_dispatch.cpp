// Full codec dispatch. Only equilibrium::io links optional file libraries.
#include "equilibrium/reader.hpp"
#include "equilibrium/writer.hpp"
#include "internal_factories.hpp"

#include <memory>

namespace equilibrium {

bool output_format_available(OutputFormat format) {
    if (format == OutputFormat::BINARY) return true;
#ifdef EQUILIBRIUM_HAVE_NETCDF
    if (format == OutputFormat::NETCDF) return true;
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
    if (format == OutputFormat::HDF5) return true;
#endif
    return false;
}

std::unique_ptr<Writer> create_writer(OutputFormat format) {
    if (format == OutputFormat::BINARY) { return make_binary_writer(); }
#ifdef EQUILIBRIUM_HAVE_NETCDF
    if (format == OutputFormat::NETCDF) { return make_netcdf_v1_writer(); }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
    if (format == OutputFormat::HDF5) { return make_hdf5_v1_writer(); }
#endif
    return nullptr;
}

std::unique_ptr<Reader> make_reader(OutputFormat format) {
    if (format == OutputFormat::BINARY) { return make_binary_reader(); }
#ifdef EQUILIBRIUM_HAVE_NETCDF
    if (format == OutputFormat::NETCDF) {
        return make_netcdf_v1_reader();  // defined in netcdf_writer.cpp
    }
#endif
#ifdef EQUILIBRIUM_HAVE_HDF5
    if (format == OutputFormat::HDF5) {
        return make_hdf5_v1_reader();  // defined in hdf5_writer.cpp
    }
#endif
    return nullptr;
}

}  // namespace equilibrium
