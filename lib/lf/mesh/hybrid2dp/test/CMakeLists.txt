include(GoogleTest)

set(sources
  mesh_factory_tests.cc
  mesh_factory_test.h
  mesh_orientation_test.cc
)

add_executable(lf.mesh.hybrid2dp.test ${sources})
target_link_libraries(lf.mesh.hybrid2dp.test PUBLIC Eigen3::Eigen Boost::boost GTest::main
  lf.mesh.hybrid2dp
  lf.mesh.test_utils
  lf.mesh.utils)
target_compile_features(lf.mesh.hybrid2dp.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.mesh.hybrid2dp.test)

# set(orisources
#   mesh_factory_test.h
#   mesh_orientation_test.cc
# )

# add_executable(lf.mesh.orientation.test ${orisources})
# target_link_libraries(lf.mesh.orientation.test PUBLIC Eigen3::Eigen Boost::boost GTest::main
#   lf.mesh.hybrid2dp
#   lf.mesh.test_utils
#   lf.mesh.utils)
# target_compile_features(lf.mesh.orientation.test PUBLIC cxx_std_17)
# gtest_discover_tests(lf.mesh.orientation.test)
