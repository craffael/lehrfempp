
include(GoogleTest)

set(sources
  gauss_jacobi_tests.cc
  gauss_legendre_tests.cc
  make_quad_rule_tests.cc
)

add_executable(lf.quad.test ${sources})
target_link_libraries(lf.quad.test PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main lf.base lf.quad)
set_target_properties(lf.quad.test PROPERTIES FOLDER "lib")
gtest_discover_tests(lf.quad.test)
