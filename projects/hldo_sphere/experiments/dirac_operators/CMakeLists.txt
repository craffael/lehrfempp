project(projects.hldo_sphere.experiments.dirac_operators)

add_library(projects.hldo_sphere.experiments.dirac_opeartor dirac_operator_experiment.h)
target_link_libraries(projects.hldo_sphere.experiments.dirac_opeartor PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh lf.mesh.utils lf.assemble lf.quad projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.experiments.dirac_opeartor PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.experiments.dirac_opeartor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(projects.hldo_sphere.experiments.dirac_opeartor PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.zero dirac_operator_zero.cc)
target_link_libraries(${PROJECT_NAME}.zero PUBLIC Eigen3::Eigen projects.hldo_sphere.experiments.dirac_opeartor projects.hldo_sphere.operators lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.zero PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.zero PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.constant dirac_operator_constant.cc)
target_link_libraries(${PROJECT_NAME}.constant PUBLIC Eigen3::Eigen projects.hldo_sphere.experiments.dirac_opeartor projects.hldo_sphere.operators lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.constant PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.constant PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.test dirac_operator_test.cc)
target_link_libraries(${PROJECT_NAME}.test PUBLIC Eigen3::Eigen projects.hldo_sphere.experiments.dirac_opeartor projects.hldo_sphere.operators lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.test PUBLIC cxx_std_17)

# copy plot file
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/plot.py ${CMAKE_CURRENT_BINARY_DIR}/plot.py COPYONLY)
