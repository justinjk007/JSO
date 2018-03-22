cmake_minimum_required (VERSION 3.2)
project (jSO)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (MSVC)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot3d.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Debug)
endif(MSVC)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/plot3d.py
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
message("-- Copied python script for plotting")

# Set compiler specific options
if (MSVC)
  add_compile_options(/MP /Gm- /W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -g)
endif(MSVC)

add_executable(run_tests tests/run_tests.cpp tests/integration_tests.cpp)
