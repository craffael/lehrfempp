
set(sources 
    curl_element_vector_provider.h
    curl_element_vector_provider.cc
    rot_whitney_one_form_dot_element_matrix_provider.h
    rot_whitney_one_form_dot_element_matrix_provider.cc
    rot_w_one_form_div_w_two_form_element_matrix_provider.h
    rot_w_one_form_div_w_two_form_element_matrix_provider.cc
    whitney_one_form_curl_element_matrix_provider.h
    whitney_one_form_curl_element_matrix_provider.cc
    whitney_one_form_grad_element_matrix_provider.h
    whitney_one_form_grad_element_matrix_provider.cc
    whitney_two_element_vector_provider.h
    whitney_two_element_vector_provider.cc
)

add_library(projects.hldo_sphere.assemble ${sources})
target_link_libraries(projects.hldo_sphere.assemble PUBLIC Eigen3::Eigen lf.base lf.geometry lf.uscalfe lf.quad lf.assemble lf.mesh lf.mesh.utils lf.io projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.assemble PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.assemble PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.hldo_sphere.assemble PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(test)
