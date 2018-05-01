
include(GoogleTest)

set(sources 
  ref_el.cc
)

add_executable(base.test ${sources})
target_link_libraries(base.test PUBLIC Eigen3::Eigen Boost::boost GTest::main base)
target_compile_features(base.test PUBLIC cxx_std_17)
gtest_discover_tests(base.test)
