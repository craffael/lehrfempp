project(projects.hldo_sphere.debugging.experiments)

add_executable(${PROJECT_NAME}.convergence_test dirac_convergence_test.cc)
target_link_libraries(${PROJECT_NAME}.convergence_test PUBLIC Eigen3::Eigen projects.hldo_sphere.debugging  lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.convergence_test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.convergence_test PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.basis_expansion_test whitney_one_basis_expansion_test.cc)
target_link_libraries(${PROJECT_NAME}.basis_expansion_test PUBLIC Eigen3::Eigen projects.hldo_sphere.debugging  lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.basis_expansion_test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.basis_expansion_test PUBLIC cxx_std_17)

add_executable(${PROJECT_NAME}.whitney_one_curl_test whitney_one_curl_convergence_test.cc)
target_link_libraries(${PROJECT_NAME}.whitney_one_curl_test PUBLIC Eigen3::Eigen projects.hldo_sphere.debugging  lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.operators)
set_target_properties(${PROJECT_NAME}.whitney_one_curl_test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.whitney_one_curl_test PUBLIC cxx_std_17)

# copy plot file
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/plotwonecurldebug.py ${CMAKE_CURRENT_BINARY_DIR}/plotwonecurldebug.py COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/plotconvdebug.py ${CMAKE_CURRENT_BINARY_DIR}/plotconvdebug.py COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/projects/hldo_sphere/post_processing/plotbasisexpdebug.py ${CMAKE_CURRENT_BINARY_DIR}/plotbasisexpdebug.py COPYONLY)
