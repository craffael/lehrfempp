
include(GoogleTest)

set(sources
  forward_iterator_tests.cc
  forward_range_tests.cc
  random_access_iterator_tests.cc
  ref_el_tests.cc
  static_vars_tests.cc
)

add_executable(lf.base.test ${sources})
target_link_libraries(lf.base.test PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.base)
target_compile_features(lf.base.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.base.test)
