// internal_factories.hpp — cross-TU declarations for the v1 NetCDF/HDF5
// backend factories. Not installed.
#ifndef EQUILIBRIUM_INTERNAL_FACTORIES_HPP_
#define EQUILIBRIUM_INTERNAL_FACTORIES_HPP_

#include "equilibrium/reader.hpp"
#include "equilibrium/writer.hpp"

#include <memory>

namespace equilibrium {

// Optional codecs are linked through equilibrium::io.
std::unique_ptr<Writer> make_netcdf_v1_writer();
std::unique_ptr<Reader> make_netcdf_v1_reader();
std::unique_ptr<Writer> make_hdf5_v1_writer();
std::unique_ptr<Reader> make_hdf5_v1_reader();

}  // namespace equilibrium

#endif  // EQUILIBRIUM_INTERNAL_FACTORIES_HPP_
