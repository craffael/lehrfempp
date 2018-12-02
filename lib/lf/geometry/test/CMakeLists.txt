
include(GoogleTest)

set(sources point_tests.cc quad_test.cc)

add_executable(lf.geometry.test ${sources})
target_link_libraries(lf.geometry.test PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.geometry)
gtest_discover_tests(lf.geometry.test)
