set(sources				
fe_space_lagrange_o1.h						
fe_tools.h
fe_tools.cc
lagr_fe.h
lagr_fe.cc
lin_fe.h
lin_fe.cc	
loc_comp_ellbvp.cc				
loc_comp_ellbvp.h
loc_comp_norms.h
loc_comp_norms.cc
mesh_function_constant.h			
mesh_function_global.h
mesh_function_traits.h
precomputed_scalar_reference_finite_element.h				
uniform_scalar_fe_space.h
uniform_scalar_fe_space.cc
uscalfe.h
)

lf_add_library(lf.uscalfe ${sources})
target_link_libraries(lf.uscalfe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad)
target_compile_features(lf.uscalfe PUBLIC cxx_std_17)

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()