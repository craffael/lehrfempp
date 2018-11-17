
include(GoogleTest)

set(sources
  gauss_jacobi_tests.cc
  gauss_legendre_tests.cc
  make_quad_rule_tests.cc
)

add_executable(quad.test ${sources})
target_link_libraries(quad.test PUBLIC Eigen3::Eigen Boost::boost GTest::main base quad)
gtest_discover_tests(quad.test)
