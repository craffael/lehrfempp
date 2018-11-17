include(GoogleTest)

set(sources
  refinement_test_utils.h
  refinement_test_utils.cc
    father_child_relation_tests.cc 
    geo_ref_test.cc
    hybrid2d_refinement_pattern_tests.cc
    regreftest.cc
)

add_executable(refinement.test ${sources})
target_link_libraries(refinement.test
  PUBLIC Eigen3::Eigen Boost::boost GTest::main
  refinement
  mesh.hybrid2d
  mesh.test_utils
  mesh.utils
  io
  io.test_utils
  )

gtest_discover_tests(refinement.test)

set(rrt_sources 
  refinement_test_utils.h
  refinement_test_utils.cc
  regreftest.cc
)

add_executable(refinement.rrt_test ${rrt_sources})
target_link_libraries(refinement.rrt_test
  PUBLIC Eigen3::Eigen Boost::boost GTest::main
  refinement
  mesh.hybrid2d
  mesh.test_utils
  mesh.utils
  io
  io.test_utils
  )

gtest_discover_tests(refinement.rrt_test)
