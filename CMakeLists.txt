cmake_minimum_required(VERSION 3.0)

project(spatial-tree)

add_compile_options(-Ofast -march=native -Wall -Wextra -Wpedantic -Werror -std=c++20 -std=gnu++20 -fverbose-asm -save-temps)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ALL_LIBRARIES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../bin/")

add_subdirectory(spatial-tree)
add_subdirectory(test)
