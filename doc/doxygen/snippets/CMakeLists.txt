# Here we want to compile all snippets to make sure that they stay up-to-date, even when the code is refactored.

set(sources
  forward_iterator.cc
  gmsh_reader.cc
  lambda_mesh_data_set.cc
  mesh_function_binary.cc
  ref_el.cc
  vtk_writer.cc
  mesh_utils.cc
)

add_library(docsnippets ${sources})
target_link_libraries(docsnippets PUBLIC Eigen3::Eigen lf.base lf.mesh)
target_compile_features(docsnippets PUBLIC cxx_std_17)
