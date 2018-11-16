set(sources
  lagr_fe.h
  lagr_fe.cc
  loc_comp_ellbvp.h
  loc_comp_ellbvp.cc
  loc_comp_norms.h
  loc_comp_norms.cc
  fe_space.h
  fe_space.cc
  fe_tools.h
  fe_tools.cc
  fe.h
)

add_library(lf.fe ${sources})
target_link_libraries(lf.fe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad)
target_compile_features(lf.fe PUBLIC cxx_std_17)

add_subdirectory(test)
