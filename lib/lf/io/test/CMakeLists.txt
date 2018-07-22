
include(GoogleTest)

set(sources
  gmsh_reader_tests.cc
  vtk_writer_tests.cc
)

add_executable(lf.io.test ${sources})
target_link_libraries(lf.io.test PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.io lf.io.test_utils lf.mesh.hybrid2d lf.mesh.hybrid2dp)
target_compile_features(lf.io.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.io.test)
