# Add a custom target to generate all the plots from the numerical experiments
add_custom_target(projects.ipdg_stokes.plots)

# Define _USE_MATH_DEFINES so that MSVC is happy
add_definitions(-D_USE_MATH_DEFINES)

# Add the different subprojects
add_subdirectory(assemble)
add_subdirectory(mesh)
add_subdirectory(post_processing)
add_subdirectory(examples)
