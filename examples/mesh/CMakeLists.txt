# cmake configuration file for compilation of LehrFEM++ examples: Mesh::Entity

set(mesh_entity_sources mesh_entity_output.cc)

add_executable(examples.mesh.mesh_entity_output ${mesh_entity_sources})


# Target link libraries. Also include lf.mesh.test? Remove lf.mesh.hybrid2d for later?

target_link_libraries(examples.mesh.mesh_entity_output PUBLIC Eigen3::Eigen
  Boost::boost GTest::main mesh.hybrid2d mesh.test_utils mesh.utils)

target_compile_features(examples.mesh.mesh_entity_output PUBLIC cxx_std_17)

set(structured_mesh_sources structured_mesh_demo.cc)
add_executable(examples.mesh.structured_mesh_demo ${structured_mesh_sources})
target_link_libraries(examples.mesh.structured_mesh_demo PUBLIC Eigen3::Eigen
  Boost::boost mesh.hybrid2d mesh.test_utils mesh.utils)

