set(sources				
fe_space_lagrange_o1.h				
fe_space_lagrange_uniform.h
fe_space_lagrange_uniform.cc			
fe_tools.h
fe_tools.cc	
lagr_fe.h	
lagr_fe.cc
lin_fe.h
lin_fe.cc
loc_comp_ellbvp.h
loc_comp_ellbvp.cc	
mesh_function_binary.h
mesh_function_constant.h
mesh_function_fe.h
mesh_function_grad_fe.h
mesh_function_global.h
mesh_function_traits.h
mesh_function_unary.h
precomputed_scalar_reference_finite_element.h				
scalar_uniform_fe_space.h
scalar_uniform_fe_space.cc
uscalfe.h
uniform_scalar_finite_element_space.h
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