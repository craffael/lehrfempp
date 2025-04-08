# Here we want to compile all snippets to make sure that they stay up-to-date, even when the code is refactored.

add_subdirectory(quick_reference)

set(sources
  assembler.cc
  coomatrix.cpp
  dofhandleruse.cc
  entityuse.cc
  fe_tools.cc
  fe_space_lagrange_o1.cc
  geometryuse.cc
  gmsh_reader.cc
  hierarchic_scalar_fe_space_snippets.cc
  line_feed_formatter.cc
  loggers.cc
  mesh_function_binary.cc
  mesh_utils.cc
  meshhierarchy.cc
  meshuse.cc
  quadruleuse.cc
  ref_el.cc
  timer.cc
  vtk_writer.cc
)

add_library(docsnippets ${sources})
target_link_libraries(docsnippets PUBLIC Eigen3::Eigen lf.base lf.mesh lf.geometry lf.assemble)
set_target_properties(docsnippets PROPERTIES FOLDER "doc")
target_compile_features(docsnippets PUBLIC cxx_std_17)
