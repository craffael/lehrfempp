# cmake configuration file for compilation of LehrFEM++ examples: Mesh::Entity

set(mesh_entity_sources mesh_entity_output.cc)

add_executable(examples.mesh.mesh_entity_output ${mesh_entity_sources})
target_link_libraries(examples.mesh.mesh_entity_output PUBLIC Eigen3::Eigen
  Boost::boost GTest::main lf.io lf.mesh.hybrid2d lf.mesh.test_utils lf.mesh.utils)
add_custom_target(examples.mesh.mesh_entity_output_run COMMAND examples.mesh.mesh_entity_output)
add_dependencies(examples_run examples.mesh.mesh_entity_output_run)

set(structured_mesh_sources structured_mesh_demo.cc)
add_executable(examples.mesh.structured_mesh_demo ${structured_mesh_sources})
target_link_libraries(examples.mesh.structured_mesh_demo PUBLIC Eigen3::Eigen
  Boost::boost lf.io lf.mesh.hybrid2d lf.mesh.test_utils lf.mesh.utils)
add_custom_target(examples.mesh.structured_mesh_demo_run COMMAND examples.mesh.structured_mesh_demo)
add_dependencies(examples_run examples.mesh.structured_mesh_demo_run)


