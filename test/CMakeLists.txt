find_package(GTest CONFIG)
if (${GTest_FOUND})
  enable_testing()
  function(create_test test_name entry_point)
    add_executable(${test_name} ${entry_point})
    target_include_directories(
      ${test_name}
      PUBLIC
      ${CMAKE_SOURCE_DIR}
     )
     target_link_libraries(${test_name} gtest gtest_main)
     add_test(${test_name} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name}")
   endfunction()
 
   create_test("test_kdtree" kdtree.cc)
   create_test("test_utils" utils.cc)
endif()