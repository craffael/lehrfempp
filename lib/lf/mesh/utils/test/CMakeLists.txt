
include(GoogleTest)

set(sources
  count_test.cc
  torus_mesh_builder_tests.cc
  tp_quad_mesh_builder_tests.cc
  tp_triag_mesh_builder_tests.cc
)

add_executable(mesh.utils.test ${sources})
target_link_libraries(mesh.utils.test
  PUBLIC Eigen3::Eigen Boost::boost GTest::main
  mesh.hybrid2d
  mesh.utils
  mesh.test_utils
)
gtest_discover_tests(mesh.utils.test)
