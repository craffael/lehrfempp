
set(sources piecewise_const_element_matrix_provider.h
	    piecewise_const_element_matrix_provider.cc
	    piecewise_const_element_vector_provider.h
	    piecewise_const_element_vector_provider.cc
	    build_system_matrix.h
	    build_system_matrix.cc
	    offset_function.h
	    offset_function.cc
)

add_library(projects.ipdg_stokes.assemble ${sources})
target_link_libraries(projects.ipdg_stokes.assemble PUBLIC Eigen3::Eigen lf.base lf.geometry lf.uscalfe lf.quad lf.assemble lf.mesh lf.mesh.utils lf.io projects.ipdg_stokes.mesh projects.ipdg_stokes.post_processing)
set_target_properties(projects.ipdg_stokes.assemble PROPERTIES FOLDER "projects/ipdg_stokes")
target_include_directories(projects.ipdg_stokes.assemble PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.ipdg_stokes.assemble PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(test)
