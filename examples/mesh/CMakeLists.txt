# cmake configuration file for compilation of LehrFEM++ examples: Mesh::Entity

set(mesh_entity_sources mesh_entity_output.cc)

add_executable(lf.examples.mesh.mesh_entity_output ${mesh_entity_sources})

# Target link libraries? See similar file for base example
