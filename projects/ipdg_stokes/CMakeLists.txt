# Add a custom target to generate all the plots from the numerical experiments
add_custom_target(projects.ipdg_stokes.plots)

# Add the different subprojects
add_subdirectory(assemble)
add_subdirectory(mesh)
add_subdirectory(post_processing)
add_subdirectory(examples)
