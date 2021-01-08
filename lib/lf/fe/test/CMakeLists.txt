include(GoogleTest)

set(src
    fe_space_hierarchic_tests.cc
)

add_executable(lf.fe.test ${src})

target_link_libraries(lf.fe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::gtest_main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.io lf.io.test_utils lf.mesh.utils lf.mesh.test_utils lf.uscalfe)
target_compile_features(lf.fe.test PUBLIC cxx_std_17)
set_target_properties(lf.fe.test PROPERTIES FOLDER "lib")
gtest_discover_tests(lf.fe.test)

