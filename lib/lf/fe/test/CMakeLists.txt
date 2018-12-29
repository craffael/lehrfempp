include(GoogleTest)

set(src
  bvp_fe_test.cc
  fe_testutils.h
  fe_testutils.cc
  full_gal_test.cc
  gfe_test.cc
  lagr_fe_test.cc
  mesh_function_traits_test.cc
  sec_ord_ell_bvp.h
  sec_ord_ell_bvp.cc
)

add_executable(lf.fe.test ${src})

target_link_libraries(lf.fe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils lf.fe)
target_compile_features(lf.fe.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.fe.test)

