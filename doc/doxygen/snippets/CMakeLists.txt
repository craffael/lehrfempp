# Here we want to compile all snippets to make sure that they stay up-to-date, even when the code is refactored.

set(sources
  assembler.cc
  fe_tools.cc
  fe_space_lagrange_o1.cc
  gmsh_reader.cc
  lambda_mesh_data_set.cc
  mesh_function_binary.cc
  mesh_utils.cc
  ref_el.cc
  vtk_writer.cc
  coomatrix.cpp
  meshhierarchy.cc
  meshuse.cc
  entityuse.cc
  geometryuse.cc
  dofhandleruse.cc
)

add_library(docsnippets ${sources})
target_link_libraries(docsnippets PUBLIC Eigen3::Eigen lf.base lf.mesh lf.geometry lf.assemble)
target_compile_features(docsnippets PUBLIC cxx_std_17)
