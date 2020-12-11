set(sources				
    fe_space_hierarchic.h
    fe_tools.h
    fe_tools.cc
    fe.h
    hierarchic_fe.h
    hierarchic_fe.cc
    loc_comp_ellbvp.h
    loc_comp_ellbvp.cc
    mesh_function_fe.h
    mesh_function_grad_fe.h
    scalar_fe_space.h
    scalar_fe_space.cc
    scalar_reference_finite_element.h
)

lf_add_library(lf.fe ${sources})
target_link_libraries(lf.fe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad)
target_compile_features(lf.fe PUBLIC cxx_std_17)

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
