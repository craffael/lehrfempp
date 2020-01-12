set(sources
  all_codim_mesh_data_set.h
  codim_mesh_data_set.h
  lambda_mesh_data_set.h
  mesh_data_set.h
  mesh_function_binary.h
  mesh_function_constant.h
  mesh_function_global.h
  mesh_function_traits.h
  mesh_function_unary.h
  print_info.cc
  print_info.h
  special_entity_sets.h
  special_entity_sets.cc
  structured_mesh_builder.h
  torus_mesh_builder.h
  torus_mesh_builder.cc
  tp_quad_mesh_builder.h
  tp_quad_mesh_builder.cc
  tp_triag_mesh_builder.h
  tp_triag_mesh_builder.cc
  utils.h
)

lf_add_library(lf.mesh.utils ${sources})
target_link_libraries(lf.mesh.utils PUBLIC Eigen3::Eigen Boost::boost lf.mesh)
target_compile_features(lf.mesh.utils PUBLIC cxx_std_17)

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
