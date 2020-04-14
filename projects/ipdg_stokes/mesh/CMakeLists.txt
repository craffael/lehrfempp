
set(sources annulus_triag_mesh_builder.h
	    annulus_triag_mesh_builder.cc
	    utils.h
	    utils.cc
)

add_library(projects.ipdg_stokes.mesh ${sources})
target_link_libraries(projects.ipdg_stokes.mesh PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh lf.mesh.utils Boost::program_options)
target_include_directories(projects.ipdg_stokes.mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.ipdg_stokes.mesh PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(test)
