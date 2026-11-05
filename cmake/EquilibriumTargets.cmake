include_guard(GLOBAL)

# Source-build entry points for alternate toolchains. Callers own compiler,
# exception, thread, and sanitizer flags; this module owns the codec sources.
function(equilibrium_add_binary_library target)
  include(${CMAKE_CURRENT_FUNCTION_LIST_DIR}/sources.cmake)
  add_library(${target} STATIC ${EQUILIBRIUM_BINARY_SOURCES})
  target_compile_features(${target} PUBLIC cxx_std_20)
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${EQUILIBRIUM_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

# Supplied targets may wrap cross-compiled libraries. Do not search the host.
function(equilibrium_add_io_library target binary_target)
  cmake_parse_arguments(PARSE_ARGV 2 codec "" "NETCDF_TARGET;HDF5_TARGET" "")
  if(codec_UNPARSED_ARGUMENTS OR codec_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR "Invalid codec target arguments")
  endif()
  get_filename_component(root "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/.." ABSOLUTE)
  add_library(${target} STATIC ${root}/src/writer_dispatch.cpp ${root}/src/io.cpp
    ${root}/src/boozer_dispatch.cpp)
  target_link_libraries(${target} PUBLIC ${binary_target})
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF)
  foreach(backend IN ITEMS NETCDF HDF5)
    if(codec_${backend}_TARGET)
      if(NOT TARGET ${codec_${backend}_TARGET})
        message(FATAL_ERROR "Missing ${backend} dependency target: ${codec_${backend}_TARGET}")
      endif()
      string(TOLOWER ${backend} source)
      target_sources(${target} PRIVATE ${root}/src/${source}_writer.cpp
        ${root}/src/boozer_${source}.cpp)
      set_property(TARGET ${target} PROPERTY EQUILIBRIUM_HAVE_${backend} ON)
      set_property(TARGET ${target} PROPERTY EQUILIBRIUM_${backend}_TARGET ${codec_${backend}_TARGET})
      target_compile_definitions(${target} PRIVATE EQUILIBRIUM_HAVE_${backend})
      target_link_libraries(${target} PRIVATE ${codec_${backend}_TARGET})
    endif()
  endforeach()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()
