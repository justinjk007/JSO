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
  add_compile_options(-Wall -Wextra -Wpedantic -g)
endif(MSVC)

# Add source code coverage option
if (CMAKE_COMPILER_IS_GNUCXX)
  option(ENABLE_COVERAGE "Enable coverage reporting for gcc" OFF)
  if (ENABLE_COVERAGE)
    message(STATUS "Enabling coverage")
    add_compile_options(--coverage -O0) # -O0 mean zero optimization
  endif()
endif()

add_executable(run_tests tests/run_tests.cpp tests/integration_tests.cpp)

# Add testing using ctest for automation
enable_testing()
add_test(testing run_tests)
