cmake_minimum_required(VERSION 3.20)
project(equilibrium_io VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EQUILIBRIUM_STANDALONE OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(EQUILIBRIUM_STANDALONE ON)
endif()
option(EQUILIBRIUM_USE_NETCDF "Read and write NetCDF containers" ON)
option(EQUILIBRIUM_USE_HDF5 "Read and write HDF5 containers" ON)
option(EQUILIBRIUM_BUILD_TESTS "Build codec and data contract tests" ${EQUILIBRIUM_STANDALONE})
option(EQUILIBRIUM_BUILD_TOOLS "Build the equilibrium converter" ${EQUILIBRIUM_STANDALONE})

include(cmake/EquilibriumTargets.cmake)
add_library(equilibrium_data INTERFACE)
add_library(equilibrium::data ALIAS equilibrium_data)
set_target_properties(equilibrium_data PROPERTIES EXPORT_NAME data)
target_compile_features(equilibrium_data INTERFACE cxx_std_20)
target_include_directories(equilibrium_data INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

equilibrium_add_binary_library(equilibrium_binary)
add_library(equilibrium::binary ALIAS equilibrium_binary)
set_target_properties(equilibrium_binary PROPERTIES EXPORT_NAME binary)
target_link_libraries(equilibrium_binary PUBLIC equilibrium::data)

set(EQUILIBRIUM_CODEC_TARGETS)
set(EQUILIBRIUM_HAVE_NETCDF OFF)
set(EQUILIBRIUM_NETCDF_PKGCONFIG OFF)
if(EQUILIBRIUM_USE_NETCDF)
  find_package(netCDF QUIET)
  if(TARGET netCDF::netCDF)
    set(EQUILIBRIUM_NETCDF_TARGET netCDF::netCDF)
  elseif(TARGET netCDF::netcdf)
    set(EQUILIBRIUM_NETCDF_TARGET netCDF::netcdf)
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(equilibrium_netcdf REQUIRED IMPORTED_TARGET netcdf)
    set(EQUILIBRIUM_NETCDF_TARGET PkgConfig::equilibrium_netcdf)
    set(EQUILIBRIUM_NETCDF_PKGCONFIG ON)
  endif()
  set(EQUILIBRIUM_HAVE_NETCDF ON)
  list(APPEND EQUILIBRIUM_CODEC_TARGETS NETCDF_TARGET ${EQUILIBRIUM_NETCDF_TARGET})
endif()

set(EQUILIBRIUM_HAVE_HDF5 OFF)
if(EQUILIBRIUM_USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C)
  set(EQUILIBRIUM_HAVE_HDF5 ON)
  list(APPEND EQUILIBRIUM_CODEC_TARGETS HDF5_TARGET HDF5::HDF5)
endif()

equilibrium_add_io_library(equilibrium_io equilibrium_binary ${EQUILIBRIUM_CODEC_TARGETS})
add_library(equilibrium::io ALIAS equilibrium_io)
set_target_properties(equilibrium_io PROPERTIES EXPORT_NAME io)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
foreach(component IN ITEMS data binary io)
  install(TARGETS equilibrium_${component} EXPORT equilibrium_${component}Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(EXPORT equilibrium_${component}Targets NAMESPACE equilibrium::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/equilibrium_io)
endforeach()
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES README.md LICENSE CHANGELOG.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/equilibrium_io)
install(FILES docs/compatibility.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/equilibrium_io/docs)
configure_package_config_file(cmake/equilibrium_ioConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/equilibrium_ioConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/equilibrium_io)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/equilibrium_ioConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/equilibrium_ioConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/equilibrium_ioConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/equilibrium_io)

if(EQUILIBRIUM_BUILD_TOOLS)
  add_executable(equilibrium-convert tools/convert.cpp)
  target_link_libraries(equilibrium-convert PRIVATE equilibrium::io)
  install(TARGETS equilibrium-convert RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(EQUILIBRIUM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
