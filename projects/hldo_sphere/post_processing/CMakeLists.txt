
set(sources 
	    norms.h
        mesh_function_whitney_zero.h
        mesh_function_whitney_one.h
        mesh_function_whitney_two.h
)

add_library(projects.hldo_sphere.post_processing ${sources})
target_link_libraries(projects.hldo_sphere.post_processing PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh lf.mesh.utils lf.assemble lf.quad projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.post_processing PROPERTIES FOLDER "projects/ipdg_stokes")
target_include_directories(projects.hldo_sphere.post_processing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.hldo_sphere.post_processing PUBLIC cxx_std_17)

# Include the tests
add_subdirectory(test)
