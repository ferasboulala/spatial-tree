cmake_minimum_required(VERSION 3.0)

project(spatial-tree)

if (USE_EMSCRIPTEN STREQUAL "1")
  set(USE_EMSCRIPTEN TRUE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -s ASYNCIFY")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -s ASYNCIFY")
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -DPLATFORM_WEB")
else()
  set(USE_EMSCRIPTEN FALSE)
endif()

# Do not use -ffast-math because it will break floating point correctness. Use it only if coordinate types are integers.
if (USE_EMSCRIPTEN)
  add_compile_options(-O3)
else()
  add_compile_options(-O3 -march=native -Wall -Wextra -fverbose-asm)
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ALL_LIBRARIES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../bin/")

if (NOT USE_EMSCRIPTEN)
  add_executable("api" examples/api.cpp)
  target_include_directories("api" PUBLIC "${CMAKE_SOURCE_DIR}/spatial-tree")
endif()

include(FetchContent)

find_package(raylib 5.0 QUIET) # QUIET or REQUIRED
if (NOT raylib_FOUND) # If there's none, fetch and build raylib
  include(FetchContent)
  FetchContent_Declare(
    raylib
    GIT_REPOSITORY https://github.com/raysan5/raylib.git
    GIT_TAG master
    GIT_PROGRESS TRUE
    GIT_SHALLOW 1
  )
  FetchContent_GetProperties(raylib)
  if (NOT raylib_POPULATED) # Have we downloaded raylib yet?
    set(FETCHCONTENT_QUIET NO)
    FetchContent_Populate(raylib)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples
    add_subdirectory(${raylib_SOURCE_DIR} ${raylib_BINARY_DIR})
  endif()
endif()

if (NOT USE_EMSCRIPTEN)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/heads/main.zip
  )
  FetchContent_Declare(
    benchmark 
    URL https://github.com/google/benchmark/archive/refs/heads/main.zip
  )
  FetchContent_MakeAvailable(googletest benchmark)

  enable_testing()
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
if (NOT USE_EMSCRIPTEN)
  create_binary("n-body" examples/n-body.cpp)
endif()
