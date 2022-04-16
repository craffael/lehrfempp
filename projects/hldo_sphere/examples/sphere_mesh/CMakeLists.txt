project(project.hldo_sphere.examples.sphere_mesh)

add_executable(${PROJECT_NAME}.main create_mesh.cc)
target_link_libraries(${PROJECT_NAME}.main PUBLIC Eigen3::Eigen lf.base lf.mesh lf.mesh.hybrid2d lf.geometry lf.io projects.hldo_sphere.mesh)
set_target_properties(${PROJECT_NAME}.main PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.main PUBLIC cxx_std_17)
