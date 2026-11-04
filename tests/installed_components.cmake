function(run)
  execute_process(COMMAND ${ARGV} RESULT_VARIABLE status
    OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(NOT status EQUAL 0)
    message(FATAL_ERROR "${ARGV}\n${output}\n${error}")
  endif()
endfunction()

set(prefix "${BINARY_DIR}/installed-test")
run("${CMAKE_COMMAND}" --install "${BINARY_DIR}" --prefix "${prefix}")
foreach(component IN ITEMS data binary io)
  set(disabled)
  if(NOT component STREQUAL "io")
    set(disabled -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON
      -DCMAKE_DISABLE_FIND_PACKAGE_netCDF=ON
      -DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON)
  endif()
  set(build "${BINARY_DIR}/consumer-${component}")
  run("${CMAKE_COMMAND}" -S "${SOURCE_DIR}/tests/consumer" -B "${build}"
    -DCOMPONENT=${component} -DCMAKE_PREFIX_PATH=${prefix}
    -DCMAKE_CXX_COMPILER=${CXX_COMPILER}
    -DCMAKE_DISABLE_FIND_PACKAGE_CUDAToolkit=ON ${disabled})
  run("${CMAKE_COMMAND}" --build "${build}")
  run("${build}/consumer")
endforeach()
