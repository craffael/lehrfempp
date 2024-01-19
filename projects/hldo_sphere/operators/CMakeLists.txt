set(sources 
    operators.h
    whitney_zero_hodge_laplacian.h
    whitney_one_hodge_laplacian.h
    whitney_two_hodge_laplacian.h
    hodge_laplacians_source_problems.h
    dirac_operator.h
    dirac_operator.cc
    dirac_operator_source_problem.h
    dirac_operator_source_problem.cc
)

add_library(projects.hldo_sphere.operators ${sources})
target_link_libraries(projects.hldo_sphere.operators PUBLIC Eigen3::Eigen lf.base lf.geometry lf.uscalfe lf.quad lf.assemble lf.mesh lf.mesh.utils lf.io projects.hldo_sphere.assemble projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.operators PROPERTIES FOLDER "projects/hldo_sphere")
target_include_directories(projects.hldo_sphere.operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(projects.hldo_sphere.operators PUBLIC cxx_std_17)


# Include the tests
add_subdirectory(test)
