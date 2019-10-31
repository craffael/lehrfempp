include(GoogleTest)

set(src 
  assembly_test_utils.h
  convection_diffusion_ell_bvp.h
  discontinuous_fe_constant_tests.cc
  discontinuous_scalar_reference_finite_element_tests.cc
  dpg_tools_tests.cc
  lagr_fe_cubic_tests.cc
  lagr_fe_quadratic_tests.cc
  lagr_test_utils.h
  loc_comp_dpg_tests.cc
  product_dof_handler_tests.cc
  trace_scalar_reference_finite_element_tests.cc)

add_executable(projects.dpg.test ${src})

target_link_libraries(projects.dpg.test PUBLIC
  Eigen3::Eigen Boost::boost GTest::gtest_main lf.mesh lf.refinement lf.assemble lf.quad
  lf.mesh.hybrid2d lf.io lf.io.test_utils lf.mesh.utils lf.mesh.test_utils lf.uscalfe projects.dpg)
target_compile_features(projects.dpg.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.dpg.test)


