
set(sources 
    debugging.h
    dirac_convergence_test.h
    dirac_convergence_test.cc
    whitney_one_basis_expansion_coeffs.h
    whitney_one_curl_test.h
)

add_library(projects.hldo_sphere.debugging ${sources})
target_link_libraries(projects.hldo_sphere.debugging PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh lf.mesh.utils lf.assemble lf.quad projects.hldo_sphere.mesh projects.hldo_sphere.operators projects.hldo_sphere.post_processing)
set_target_properties(projects.hldo_sphere.debugging PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.debugging PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(projects.hldo_sphere.debugging PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(experiments)
