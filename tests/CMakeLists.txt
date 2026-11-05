option(EQUILIBRIUM_SANITIZER_TESTS "Build ASan/UBSan codec tests" OFF)
option(EQUILIBRIUM_PYTHON_TESTS "Test the Python reader (requires NumPy and enabled codec extras)" OFF)
option(EQUILIBRIUM_JS_TESTS "Test JavaScript payload inspection (requires Node)" OFF)

# Tests that inspect or corrupt native containers need their C APIs too.
function(equilibrium_test_backends target)
  if(EQUILIBRIUM_HAVE_NETCDF)
    target_compile_definitions(${target} PRIVATE EQUILIBRIUM_HAVE_NETCDF)
    target_link_libraries(${target} PRIVATE ${EQUILIBRIUM_NETCDF_TARGET})
  endif()
  if(EQUILIBRIUM_HAVE_HDF5)
    target_compile_definitions(${target} PRIVATE EQUILIBRIUM_HAVE_HDF5)
    target_link_libraries(${target} PRIVATE ${EQUILIBRIUM_HDF5_TARGET})
  endif()
endfunction()

function(equilibrium_sanitize target)
  target_compile_options(${target} PRIVATE
    -fsanitize=address,undefined -fno-omit-frame-pointer
    -fno-sanitize-recover=undefined)
  target_link_options(${target} PRIVATE -fsanitize=address,undefined)
endfunction()

function(equilibrium_add_test name source library sanitized)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE ${library})
  equilibrium_test_backends(${name})
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  endif()
  set(labels "unit;io")
  if(sanitized)
    equilibrium_sanitize(${name})
    list(APPEND labels sanitizer asan)
  endif()
  # Isolate fixtures even when normal and sanitized tests run in parallel.
  set(directory ${CMAKE_CURRENT_BINARY_DIR}/${name}-data)
  file(MAKE_DIRECTORY ${directory})
  add_test(NAME ${name} COMMAND ${name} ${directory})
  set_tests_properties(${name} PROPERTIES
    WORKING_DIRECTORY ${directory} LABELS "${labels}" TIMEOUT 120)
endfunction()

set(codec_tests
  host_io checkpoint asymmetric_io io_golden output_failure spectral_payload
  boozer_io)
if(EQUILIBRIUM_HAVE_NETCDF OR EQUILIBRIUM_HAVE_HDF5)
  list(APPEND codec_tests io_restart_offsets)
  # Sparse-file fixtures with 2^40 extents require a native 64-bit filesystem.
  # The standalone converter tests cover malformed Wasm inputs separately.
  if(NOT EMSCRIPTEN)
    list(APPEND codec_tests io_malformed_shapes)
  endif()
endif()

equilibrium_add_test(equilibrium_contract contract.cpp equilibrium::io FALSE)
foreach(test IN LISTS codec_tests)
  equilibrium_add_test(equilibrium_${test} test_${test}.cpp equilibrium::io FALSE)
endforeach()

if(EQUILIBRIUM_SANITIZER_TESTS)
  if(EMSCRIPTEN OR NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    message(FATAL_ERROR "Sanitizer tests require a native GCC or Clang build")
  endif()
  equilibrium_add_binary_library(equilibrium_binary_asan)
  equilibrium_add_io_library(equilibrium_io_asan equilibrium_binary_asan ${EQUILIBRIUM_CODEC_TARGETS})
  equilibrium_sanitize(equilibrium_binary_asan)
  equilibrium_sanitize(equilibrium_io_asan)
  equilibrium_add_test(equilibrium_contract_asan contract.cpp equilibrium_io_asan TRUE)
  foreach(test IN LISTS codec_tests)
    equilibrium_add_test(equilibrium_${test}_asan test_${test}.cpp equilibrium_io_asan TRUE)
  endforeach()
endif()

if(EQUILIBRIUM_STANDALONE AND NOT CMAKE_CROSSCOMPILING)
  add_test(NAME equilibrium_installed_components
    COMMAND ${CMAKE_COMMAND}
      -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
      -DBINARY_DIR=${PROJECT_BINARY_DIR}
      -DCXX_COMPILER=${CMAKE_CXX_COMPILER}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/installed_components.cmake)
endif()

if(EQUILIBRIUM_PYTHON_TESTS)
  if(CMAKE_CROSSCOMPILING)
    message(FATAL_ERROR "Python fixtures require a native build")
  endif()
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME equilibrium_python_boozer
    COMMAND ${CMAKE_COMMAND} -E env
      PYTHONPATH=${PROJECT_SOURCE_DIR}/python
      PYTHONPYCACHEPREFIX=${CMAKE_CURRENT_BINARY_DIR}/pycache
      TMPDIR=${CMAKE_CURRENT_BINARY_DIR}
      ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/test_boozer_io.py)
  add_test(NAME equilibrium_python
    COMMAND ${CMAKE_COMMAND} -E env
      PYTHONPATH=${PROJECT_SOURCE_DIR}/python
      PYTHONPYCACHEPREFIX=${CMAKE_CURRENT_BINARY_DIR}/pycache
      TMPDIR=${CMAKE_CURRENT_BINARY_DIR}
      ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/test_state_io.py
      $<TARGET_FILE:equilibrium_asymmetric_io>)
endif()

if(EQUILIBRIUM_JS_TESTS)
  find_program(EQUILIBRIUM_NODE NAMES node nodejs REQUIRED)
  add_test(NAME equilibrium_js
    COMMAND ${EQUILIBRIUM_NODE} ${CMAKE_CURRENT_SOURCE_DIR}/test_binary_payload.mjs)
endif()
