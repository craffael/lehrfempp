# Add a custom target to generate all the plots from the numerical experiments
add_custom_target(projects.hldo_sphere.plots)
set_target_properties(projects.hldo_sphere.plots PROPERTIES FOLDER "projects/hldo_sphere")

# Add the different subprojects
add_subdirectory(mesh)
add_subdirectory(examples)
add_subdirectory(assemble)
