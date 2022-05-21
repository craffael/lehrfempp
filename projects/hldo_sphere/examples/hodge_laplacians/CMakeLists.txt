project(projects.hldo_sphere.examples.hodge_laplacians)

add_executable(${PROJECT_NAME}.main hodge_laplacians.cc)
target_link_libraries(${PROJECT_NAME}.main PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.main PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.main PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.zero hodge_laplacian_zero.cc)
target_link_libraries(${PROJECT_NAME}.zero PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.zero PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.zero PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.constant hodge_laplacian_constant.cc)
target_link_libraries(${PROJECT_NAME}.constant PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.constant PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.constant PUBLIC cxx_std_17)

# copy plot file
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/plot.py ${CMAKE_CURRENT_BINARY_DIR}/plot.py COPYONLY)
