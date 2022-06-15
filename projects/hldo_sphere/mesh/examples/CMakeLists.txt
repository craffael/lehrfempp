project(project.hldo_sphere.mesh.examples)

add_executable(${PROJECT_NAME}.main create_mesh.cc)
target_link_libraries(${PROJECT_NAME}.main PUBLIC Eigen3::Eigen lf.base lf.mesh lf.mesh.hybrid2d lf.geometry lf.io projects.hldo_sphere.mesh)
set_target_properties(${PROJECT_NAME}.main PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.main PUBLIC cxx_std_17)

# create direcotry for the meshes
add_custom_target(build-time-make-directory ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory meshs)
