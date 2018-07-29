include(GoogleTest)

set(sources refine_tests.cc geo_ref_test.cc regreftest.cc)

add_executable(lf.refinement.test ${sources})
target_link_libraries(lf.refinement.test
  PUBLIC Eigen3::Eigen Boost::boost GTest::main
  lf.refinement
  lf.mesh.hybrid2dp
  lf.mesh.test_utils
  lf.mesh.utils
  lf.io)
target_compile_features(lf.refinement.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.refinement.test)
