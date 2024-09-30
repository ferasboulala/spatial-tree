cmake_minimum_required(VERSION 3.0)

project(spatial-tree)

# Set compiler and linker flags for web.
if (USE_EMSCRIPTEN STREQUAL "1")
  set(USE_EMSCRIPTEN TRUE)
else()
  set(USE_EMSCRIPTEN FALSE)
endif()

# Set compiler flags.
if (USE_EMSCRIPTEN)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -s ASYNCIFY")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -s ASYNCIFY")
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -DPLATFORM_WEB")
  add_compile_options(-O3)
else()
  # Do not use -ffast-math because it will break floating point correctness. Use it only if coordinate types are integers.
  add_compile_options(-O3 -march=native -Wall -Wextra -fverbose-asm)
endif()

# Set C++ standard and output directories.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ALL_LIBRARIES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../bin/")

# Build the API example when not compiling for web.
if (NOT USE_EMSCRIPTEN)
  add_executable("api" examples/api.cpp)
  target_include_directories("api" PUBLIC "${CMAKE_SOURCE_DIR}/spatial-tree")
endif()

set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
add_subdirectory(third-party/raylib)

if (NOT USE_EMSCRIPTEN)
  enable_testing()
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
  add_subdirectory(third-party/benchmark)
  add_subdirectory(third-party/googletest)

  function(create_test test_name entry_point)
    add_executable(${test_name} ${entry_point})
    target_include_directories(${test_name} PUBLIC ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    target_include_directories(
      ${test_name}
      PUBLIC
      ${CMAKE_SOURCE_DIR}
      # Temporary, for Mac.
      /usr/local/include/
     )
     target_link_libraries(${test_name} gtest gtest_main)
     add_test(${test_name} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name}")
   endfunction()

   create_test("utest" spatial-tree/utest.cpp)

  function(create_benchmark benchmark_name entry_point)
    add_executable(${benchmark_name} ${entry_point})
    target_include_directories(
      ${benchmark_name}
      PUBLIC
      ${CMAKE_SOURCE_DIR}
      ${benchmark_SOURCE_DIR}/include
     )
     target_link_libraries(${benchmark_name} benchmark)
   endfunction()
 
   create_benchmark("bench" spatial-tree/bench.cpp)
endif()

function(create_binary target_name entry_point)
  add_executable(${target_name} ${entry_point})
  target_include_directories(
    ${target_name}
    PUBLIC
    ${CMAKE_SOURCE_DIR}/spatial-tree)
  target_link_libraries(
    ${target_name}
    raylib
    )
endfunction()

create_binary("draw" examples/draw.cpp)
