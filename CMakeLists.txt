cmake_minimum_required (VERSION 3.2)
project (jSO)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set C++11 as the required standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler specific options
if (MSVC)
  add_compile_options(/MP /Gm- /W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -g -O0)
endif(MSVC)

add_executable(run_tests tests/run_tests.cpp tests/integration_tests.cpp)

# Add source code coverage option
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_BUILD_TYPE Debug)
  option(ENABLE_COVERAGE "Enable coverage reporting for gcc/clang" OFF)
  if (ENABLE_COVERAGE)
    message(STATUS "-- Enabling coverage")
    include(CodeCoverage)
    setup_target_for_coverage(NAME ${PROJECT_NAME}_coverage EXECUTABLE run_tests -j 4 DEPENDENCIES run_tests)
  endif()
endif()
