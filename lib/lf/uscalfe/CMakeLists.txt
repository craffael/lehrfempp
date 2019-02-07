set(sources
loc_comp_ellbvp.cc				
loc_comp_ellbvp.h
fe_space.h					
loc_comp_norms.cc
fe_space_lagrange_o1.h				
loc_comp_norms.h
fe_space_lagrange_uniform.cc			
mesh_function_constant.h
fe_space_lagrange_uniform.h			
mesh_function_global.h
fe_tools.cc					
mesh_function_traits.h
fe_tools.h					precomputed_scalar_reference_finite_element.h
lagr_fe.cc					
scalar_uniform_fe_space.cc
lagr_fe.h					
scalar_uniform_fe_space.h
lagrfe.h					
lin_fe.cc					
uniform_scalar_finite_element_space.h
lin_fe.h					
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