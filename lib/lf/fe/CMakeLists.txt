set(sources				
    fe_point.h
    fe_tools.h
    fe_tools.cc
    fe.h
    hierarchic_fe.h
    hierarchic_fe.cc
    hierarchic_scalar_fe_space.h
    loc_comp_ellbvp.h
    loc_comp_ellbvp.cc
    mesh_function_fe.h
    mesh_function_grad_fe.h
    prolongation.h
    scalar_fe_space.h
    scalar_fe_space.cc
    scalar_reference_finite_element.h
)

lf_add_library(lf.fe ${sources})
target_link_libraries(lf.fe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble lf.quad)
target_compile_features(lf.fe PUBLIC cxx_std_17)
set_target_properties(lf.fe PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()

add_subdirectory(test_utils)