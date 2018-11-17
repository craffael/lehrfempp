
include(GoogleTest)

set(sources
  gmsh_reader_tests.cc
  vtk_writer_tests.cc
)

add_executable(io.test ${sources})
target_link_libraries(io.test PUBLIC Eigen3::Eigen Boost::boost GTest::main io io.test_utils mesh.hybrid2d mesh.test_utils)
gtest_discover_tests(io.test)
