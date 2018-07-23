# cmake configuration file for compilation of LehrFEM++ examples: Mesh::Entity

set(mesh_entity_sources mesh_entity_output.cc)

add_executable(lf.examples.mesh.mesh_entity_output ${mesh_entity_sources})

# Target link libraries. Also include lf.mesh.test? Remove lf.mesh.hybrid2d for later?
target_link_libraries(lf.examples.mesh.mesh_entity_output PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.hybrid2d lf.mesh.hybrid2dp lf.mesh.test_utils lf.mesh.utils)

target_compile_features(lf.examples.mesh.mesh_entity_output PUBLIC cxx_std_17)
