project(projects.hldo_sphere.examples.dirac_operators)

add_executable(${PROJECT_NAME}.zero dirac_operator_zero.cc)
target_link_libraries(${PROJECT_NAME}.zero PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.zero PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.zero PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.constant dirac_operator_constant.cc)
target_link_libraries(${PROJECT_NAME}.constant PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.constant PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.constant PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.test dirac_operator_test.cc)
target_link_libraries(${PROJECT_NAME}.test PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.test PUBLIC cxx_std_17)
