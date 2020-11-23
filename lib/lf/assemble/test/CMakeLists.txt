include(GoogleTest)

set(sources
  assembly_tests.cc
  coomatrix_tests.cc
)

add_executable(lf.assemble.test ${sources})
target_link_libraries(lf.assemble.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::gtest_main lf.mesh lf.assemble
  lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils)
target_compile_features(lf.assemble.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.assemble.test)

set_target_properties(lf.assemble.test PROPERTIES FOLDER "lib")
