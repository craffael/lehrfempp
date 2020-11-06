
include(GoogleTest)

set(sources
  eigen_tools_tests.cc
  ref_el_tests.cc
  timer_tests.cc
)

add_executable(lf.base.test ${sources})
target_link_libraries(lf.base.test PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main lf.base)
set_target_properties(lf.base.test PROPERTIES FOLDER "lib")
gtest_discover_tests(lf.base.test)
