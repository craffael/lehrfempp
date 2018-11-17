
include(GoogleTest)

set(sources point_tests.cc quad_test.cc)

add_executable(geometry.test ${sources})
target_link_libraries(geometry.test PUBLIC Eigen3::Eigen Boost::boost GTest::main geometry)
gtest_discover_tests(geometry.test)
