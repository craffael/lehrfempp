set(sources
  fe_space_lagrange_o1.h
  fe_space_lagrange_uniform.cc
  fe_space_lagrange_uniform.h
  fe_tools.cc
  fe_tools.h
  lagrfe.h
  mesh_function_traits.h
  lagr_fe.h
  lagr_fe.cc
  lin_fe.h
  lin_fe.cc
  loc_comp_ellbvp.h
  loc_comp_ellbvp.cc
  loc_comp_norms.h
  loc_comp_norms.cc
  precomputed_scalar_reference_finite_element.h
)

add_library(lf.lagrfe ${sources})
target_link_libraries(lf.lagrfe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad)
target_compile_features(lf.lagrfe PUBLIC cxx_std_17)

add_subdirectory(test)
