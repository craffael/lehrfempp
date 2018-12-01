include(GoogleTest)

add_executable(lf.fe.test lagr_fe_test.cc)

target_link_libraries(lf.fe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils lf.fe)
target_compile_features(lf.fe.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.fe.test)

add_executable(lf.fe.gfetest gfe_test.cc)
target_link_libraries(lf.fe.gfetest PUBLIC
  Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils lf.fe)
target_compile_features(lf.fe.gfetest PUBLIC cxx_std_17)
gtest_discover_tests(lf.fe.gfetest)

add_executable(lf.fe.fullgaltest full_gal_test.cc)
target_link_libraries(lf.fe.fullgaltest PUBLIC
  Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils lf.fe)
target_compile_features(lf.fe.fullgaltest PUBLIC cxx_std_17)
gtest_discover_tests(lf.fe.fullgaltest)
