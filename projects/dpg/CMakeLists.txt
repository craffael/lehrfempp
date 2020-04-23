set(src 
  discontinuous_fe_constant.h
  discontinuous_scalar_reference_finite_element.h
  dpg.h
  dpg_element_matrix_provider.h
  dpg_element_vector_provider.h
  dpg_tools.cc
  dpg_tools.h
  loc_comp_dpg.h
  product_dofhandler.cc
  product_dofhandler.h
  product_element_matrix_provider.h
  product_element_matrix_provider_builder.h
  product_element_vector_provider.h
  product_element_vector_provider_builder.h
  product_fe_space.h
  product_fe_space_factory.h
  sub_element_matrix_provider.h
  sub_element_vector_provider.h
  trace_scalar_reference_finite_element.h)

add_library(projects.dpg ${src})
target_link_libraries(projects.dpg  
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.quad
  lf.assemble
  lf.io
  lf.uscalfe
  )

target_include_directories(projects.dpg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(projects.dpg PUBLIC cxx_std_17)


#add tests
if (LF_ENABLE_TESTING)
  add_subdirectory(test) 
endif()

#add convergence studies
add_subdirectory(convergence_studies)
