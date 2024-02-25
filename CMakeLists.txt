cmake_minimum_required(VERSION 3.0)

project(spatial-tree)

# Do not use -ffast-math because it will break floating point correctness. Use it only if coordinate types are integers.
add_compile_options(-O3 -march=native -Wall -Wextra -fverbose-asm -save-temps)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(ALL_LIBRARIES)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../bin/")

find_package(GTest CONFIG)
if (${GTest_FOUND})
  enable_testing()
  function(create_test test_name entry_point)
    add_executable(${test_name} ${entry_point})
    target_include_directories(
      ${test_name}
      PUBLIC
      ${CMAKE_SOURCE_DIR}
      # Temporary, for Mac.
      /usr/local/include/
     )
     target_link_libraries(${test_name} GTest::gtest_main)
     add_test(${test_name} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name}")
   endfunction()

   create_test("utest" src/utest.cpp)
endif()

find_package(benchmark CONFIG)
if (${benchmark_FOUND})
  function(create_benchmark benchmark_name entry_point)
    add_executable(${benchmark_name} ${entry_point})
    target_include_directories(
      ${benchmark_name}
      PUBLIC
      ${CMAKE_SOURCE_DIR}
     )
     target_link_libraries(${benchmark_name} benchmark::benchmark)
   endfunction()
 
   create_benchmark("bench" src/bench.cpp)
endif()

add_executable("example" src/example.cpp)

find_package(OpenCV)
if (DEFINED OpenCV_LIBRARIES)
  function(create_binary target_name entry_point)
    add_executable(${target_name} ${entry_point})
    target_include_directories(
      ${target_name}
      PUBLIC
      ${CMAKE_SOURCE_DIR}
      ${OpenCV_INCLUDE_DIRS})
    target_link_libraries(
      ${target_name}
      ${OpenCV_LIBRARIES})
  endfunction()
  create_binary("interactive" src/interactive.cpp)
endif()
