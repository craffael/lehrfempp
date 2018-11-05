include(GoogleTest)

set(sources
  lagr_fe_test.cc
)

add_executable(lf.fe.test ${sources})
target_link_libraries(lf.fe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils lf.fe)
target_compile_features(lf.fe.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.fe.test)
