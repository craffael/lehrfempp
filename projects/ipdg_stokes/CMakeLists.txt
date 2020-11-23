# Add a custom target to generate all the plots from the numerical experiments
add_custom_target(projects.ipdg_stokes.plots)
set_target_properties(projects.ipdg_stokes.plots PROPERTIES FOLDER "projects/ipdg_stokes")

# Add the different subprojects
add_subdirectory(assemble)
add_subdirectory(mesh)
add_subdirectory(post_processing)
add_subdirectory(examples)
