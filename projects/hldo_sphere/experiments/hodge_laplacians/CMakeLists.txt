project(projects.hldo_sphere.experiments.hodge_laplacians)

add_library(projects.hldo_sphere.experiments.hodge_laplacian hodge_laplacian_experiment.h)
target_link_libraries(projects.hldo_sphere.experiments.hodge_laplacian PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh lf.mesh.utils lf.assemble lf.quad projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.experiments.hodge_laplacian PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.experiments.hodge_laplacian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(projects.hldo_sphere.experiments.hodge_laplacian PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.zero hodge_laplacian_zero.cc)
target_link_libraries(${PROJECT_NAME}.zero PUBLIC Eigen3::Eigen projects.hldo_sphere.experiments.hodge_laplacian projects.hldo_sphere.operators lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.zero PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.zero PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.test hodge_laplacian_test.cc)
target_link_libraries(${PROJECT_NAME}.test PUBLIC Eigen3::Eigen projects.hldo_sphere.experiments.hodge_laplacian projects.hldo_sphere.operators lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.test PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.constant hodge_laplacian_constant.cc)
target_link_libraries(${PROJECT_NAME}.constant PUBLIC Eigen3::Eigen projects.hldo_sphere.experiments.hodge_laplacian projects.hldo_sphere.operators lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.constant PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.constant PUBLIC cxx_std_17)

# copy plot file
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/plotlaplace.py ${CMAKE_CURRENT_BINARY_DIR}/plotlaplace.py COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/ploteigenvalues.py ${CMAKE_CURRENT_BINARY_DIR}/ploteigenvalues.py COPYONLY)
