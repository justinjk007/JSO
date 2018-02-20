cmake_minimum_required (VERSION 3.2)
project (JSO)

# Set compiler specific options
if (MSVC)
  add_compile_options(/MP /Gm- /W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic -g)
endif(MSVC)

# add_executable(main main.cc cec17_test_func.cc lshade.cc search_algorithm.cc)
