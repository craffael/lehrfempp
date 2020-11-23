set(sources				
fe_space_lagrange_o1.h
fe_space_lagrange_o2.h
fe_space_lagrange_o3.h						
lagr_fe.h
lagr_fe.cc
lin_fe.h
lin_fe.cc
loc_comp_ellbvp.h
loc_comp_ellbvp.cc
precomputed_scalar_reference_finite_element.h				
uniform_scalar_fe_space.h
uniform_scalar_fe_space.cc
uscalfe.h
uscalfe.h
prolongation.h
)

lf_add_library(lf.uscalfe ${sources})
target_link_libraries(lf.uscalfe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad lf.fe)
set_target_properties(lf.uscalfe PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
