include(GoogleTest)

set(sources
  assembly_tests.cc
)

add_executable(assemble.test ${sources})
target_link_libraries(assemble.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::main mesh assemble
  mesh.hybrid2d mesh.utils mesh.test_utils)
target_compile_features(assemble.test PUBLIC cxx_std_17)
gtest_discover_tests(assemble.test)
