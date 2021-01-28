include(GoogleTest)

set(src
  fe_testutils.h
  fe_testutils.cc
  sec_ord_ell_bvp.h
  #sec_ord_ell_bvp.cc // nothing in the cc file, all is in the header
  loc_comp_test.cc
  bvp_fe_tests.cc
  full_gal_tests.cc
  gfe_tests.cc
  lagr_fe_tests.cc
)

add_executable(lf.uscalfe.test ${src})

target_link_libraries(lf.uscalfe.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::gtest_main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.io lf.io.test_utils lf.mesh.utils lf.mesh.test_utils lf.uscalfe)
set_target_properties(lf.uscalfe.test PROPERTIES FOLDER "lib")
gtest_discover_tests(lf.uscalfe.test)

if(WIN32)
    target_compile_options(lf.uscalfe.test PRIVATE /bigobj)
endif()

