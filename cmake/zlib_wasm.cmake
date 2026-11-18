# Emscripten emits static archives for both zlib 1.3.1 targets. Give the
# unused shared target a distinct directory so Ninja can generate its rules.
cmake_language(DEFER CALL set_target_properties zlib PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/shared")
