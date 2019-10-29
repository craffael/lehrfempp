set(sources
  check_child_geometry.cc
  check_child_geometry.h
  check_integration_element.cc
  check_integration_element.h
  check_jacobian.cc
  check_jacobian.h
  check_jacobian_inverse_gramian.cc
  check_jacobian_inverse_gramian.h
  check_sub_geometry.cc
  check_sub_geometry.h
  test_utils.cc
  test_utils.h
)

lf_add_library(lf.geometry.test_utils ${sources})
target_link_libraries(lf.geometry.test_utils PUBLIC
  GTest::gtest_main
  lf.base
  lf.geometry
  lf.quad
  lf.refinement
)