
include(GoogleTest)

set(sources
  mesh_tests.cc
)

add_executable(lf.mesh.hybrid2d.test ${sources})
target_link_libraries(lf.mesh.hybrid2d.test PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.hybrid2d lf.mesh.test_utils)
target_compile_features(lf.mesh.hybrid2d.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.mesh.hybrid2d.test)
