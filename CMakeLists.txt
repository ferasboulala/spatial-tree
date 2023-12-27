cmake_minimum_required(VERSION 3.0)

project(spatial-tree)

# Do not use -ffast-math because it will break floating point assumption.
add_compile_options(-O3 -march=native -Wall -Wextra -Wpedantic -Werror -std=c++20 -std=gnu++20 -fverbose-asm -save-temps)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ALL_LIBRARIES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../bin/")

add_subdirectory(spatial-tree)
add_subdirectory(test)
add_subdirectory(benchmark)
