# Each toolchain variant compiles this same source list with its own flags.
get_filename_component(EQUILIBRIUM_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(EQUILIBRIUM_BINARY_SOURCES
  ${EQUILIBRIUM_SOURCE_DIR}/src/boozer_binary.cpp
  ${EQUILIBRIUM_SOURCE_DIR}/src/boozer_common.cpp
  ${EQUILIBRIUM_SOURCE_DIR}/src/spectral_payload.cpp
  ${EQUILIBRIUM_SOURCE_DIR}/src/binary_io.cpp
  ${EQUILIBRIUM_SOURCE_DIR}/src/output_spec.cpp
  ${EQUILIBRIUM_SOURCE_DIR}/src/versioned_binary.cpp
  ${EQUILIBRIUM_SOURCE_DIR}/src/checkpoint.cpp)
