include(GoogleTest)

set(sources
  refinement_test_utils.h
  refinement_test_utils.cc
    father_child_relation_tests.cc 
    geo_ref_test.cc
    hybrid2d_refinement_pattern_tests.cc
    regreftest.cc
    mesh_function_transfer_tests.cc
)

add_executable(lf.refinement.test ${sources})
target_link_libraries(lf.refinement.test
  PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main
  lf.refinement
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.io
  lf.io.test_utils
  )
set_target_properties(lf.refinement.test PROPERTIES FOLDER "lib")

gtest_discover_tests(lf.refinement.test)

set(rrt_sources 
  refinement_test_utils.h
  refinement_test_utils.cc
  regreftest.cc
)

add_executable(lf.refinement.rrt_test ${rrt_sources})
target_link_libraries(lf.refinement.rrt_test
  PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main
  lf.refinement
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.io
  lf.io.test_utils
  )
set_target_properties(lf.refinement.rrt_test PROPERTIES FOLDER "lib")

gtest_discover_tests(lf.refinement.rrt_test)
