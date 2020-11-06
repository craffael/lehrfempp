
set(sources solution_to_mesh_data_set.h
	    solution_to_mesh_data_set.cc
	    norms.h
)

add_library(projects.ipdg_stokes.post_processing ${sources})
target_link_libraries(projects.ipdg_stokes.post_processing PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh lf.mesh.utils lf.assemble lf.quad projects.ipdg_stokes.mesh)
set_target_properties(projects.ipdg_stokes.post_processing PROPERTIES FOLDER "projects/ipdg_stokes")
target_include_directories(projects.ipdg_stokes.post_processing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.ipdg_stokes.post_processing PUBLIC cxx_std_17)

# Include the tests
add_subdirectory(test)
