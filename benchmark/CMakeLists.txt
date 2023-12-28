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
 
   create_benchmark("benchmark_quadtree" quadtree.cc)
   create_benchmark("benchmark_kdtree" kdtree.cc)
   create_benchmark("benchmark_hashtable" hashtable.cc)
endif()