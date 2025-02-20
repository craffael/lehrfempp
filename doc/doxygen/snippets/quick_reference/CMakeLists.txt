# Here we want to compile all snippets to make sure that they stay up-to-date, even when the code is refactored.

set(sources
    qr_fe_space_snippets.cpp
    qr_geometry_snippets.cpp
    qr_quad_snippets.cpp
    qr_mesh_snippets.cpp
    qr_bc_snippets.cpp
    qr_dofh_snippets.cpp
    
)

add_library(docsnippets_qr ${sources})
target_link_libraries(docsnippets_qr PUBLIC Eigen3::Eigen lf.base lf.mesh lf.geometry lf.assemble)
set_target_properties(docsnippets_qr PROPERTIES FOLDER "doc")
target_compile_features(docsnippets_qr PUBLIC cxx_std_17)
