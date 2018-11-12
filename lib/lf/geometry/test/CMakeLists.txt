
include(GoogleTest)

set(sources geometry_tests.cc point_tests.cc quad_test.cc)

add_executable(lf.geometry.test ${sources})
target_link_libraries(lf.geometry.test
  PUBLIC Eigen3::Eigen Boost::boost GTest::main
  lf.geometry
  lf.quad
)

target_compile_features(lf.geometry.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.geometry.test)
