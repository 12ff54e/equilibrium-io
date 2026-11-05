# netCDF 4.9.3 bundles an older FindZLIB that does not create ZLIB::ZLIB.
# HDF5's exported static target needs it during netCDF's configure-time link
# probes. Use CMake's module with the explicitly pinned Wasm library/header.
if(NOT TARGET ZLIB::ZLIB)
  include("${CMAKE_ROOT}/Modules/FindZLIB.cmake")
endif()
if(NOT TARGET ZLIB::ZLIB)
  message(FATAL_ERROR "The pinned Wasm zlib library was not found")
endif()
