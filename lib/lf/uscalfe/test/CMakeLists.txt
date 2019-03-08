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
  loc_comp_test.cc
)

add_executable(lf.uscalfe.test ${src})

target_link_libraries(lf.uscalfe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.mesh.utils lf.mesh.test_utils lf.uscalfe)
target_compile_features(lf.uscalfe.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.uscalfe.test)

