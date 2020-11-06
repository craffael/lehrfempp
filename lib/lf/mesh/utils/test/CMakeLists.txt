
include(GoogleTest)

set(sources
  count_test.cc
  mesh_function_utils.h
  mesh_function_traits_tests.cc
  mesh_function_binary_tests.cc
  mesh_function_unary_tests.cc
  torus_mesh_builder_tests.cc
  tp_quad_mesh_builder_tests.cc
  tp_triag_mesh_builder_tests.cc
)

add_executable(lf.mesh.utils.test ${sources})
target_link_libraries(lf.mesh.utils.test
  PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main
  lf.io
  lf.io.test_utils
  lf.mesh.hybrid2d
  lf.mesh.utils
  lf.mesh.test_utils
  lf.quad
)
set_target_properties(lf.mesh.utils.test PROPERTIES FOLDER "lib")
gtest_discover_tests(lf.mesh.utils.test)
