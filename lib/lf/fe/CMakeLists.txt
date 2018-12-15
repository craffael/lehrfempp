set(sources
  fe_space.cc
  fe_space.h
  fe_testutils.h
  fe_testutils.cc
  fe_tools.cc
  fe_tools.h
  fe.h
  is_mesh_function.h
  lagr_fe.h
  lagr_fe.cc
  lin_fe.h
  lin_fe.cc
  loc_comp_ellbvp.h
  loc_comp_ellbvp.cc
  loc_comp_norms.h
  loc_comp_norms.cc
  sec_ord_ell_bvp.h
  sec_ord_ell_bvp.cc
)

add_library(lf.fe ${sources})
target_link_libraries(lf.fe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad)
target_compile_features(lf.fe PUBLIC cxx_std_17)

add_subdirectory(test)
