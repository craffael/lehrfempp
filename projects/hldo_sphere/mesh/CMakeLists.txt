
set(sources sphere_triag_mesh_builder.h
	    sphere_triag_mesh_builder.cc
)

add_library(projects.hldo_sphere.mesh ${sources})
target_link_libraries(projects.hldo_sphere.mesh PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh.hybrid2d lf.mesh lf.mesh.utils Boost::program_options)
set_target_properties(projects.hldo_sphere.mesh PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.hldo_sphere.mesh PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(test)

