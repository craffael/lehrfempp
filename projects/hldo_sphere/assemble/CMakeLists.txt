
set(sources 
    laplace_matrix_provider.h
    laplace_matrix_provider.cc
    mass_matrix_provider.h
    mass_matrix_provider.cc
    load_vector_provider.h
    whitney_one_curl_curl_matrix_provider.h
    whitney_one_curl_curl_matrix_provider.cc
    whitney_one_grad_matrix_provider.h
    whitney_one_grad_matrix_provider.cc
    whitney_one_mass_matrix_provider.h
    whitney_one_mass_matrix_provider.cc
    whitney_one_vector_provider.h
    whitney_two_mass_matrix_provider.h
    whitney_two_mass_matrix_provider.cc
    rot_whitney_one_div_matrix_provider.h
    rot_whitney_one_div_matrix_provider.cc
    whitney_two_vector_provider.h
    whitney_two_vector_provider.cc
)

add_library(projects.hldo_sphere.assemble ${sources})
target_link_libraries(projects.hldo_sphere.assemble PUBLIC Eigen3::Eigen lf.base lf.geometry lf.uscalfe lf.quad lf.assemble lf.mesh lf.mesh.utils lf.io projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.assemble PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.assemble PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.hldo_sphere.assemble PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(test)
