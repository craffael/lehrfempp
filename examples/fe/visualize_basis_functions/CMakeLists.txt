# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly
set(target examples.fe.visualize_basis_functions)
add_executable(${target} visualize_basis_functions.cc)
target_link_libraries(${target}
  PUBLIC
  lf.fe
  lf.io
  lf.mesh.hybrid2d
  )
set_target_properties(${target} PROPERTIES FOLDER "examples")
add_custom_target("${target}_run" COMMAND ${target})
add_dependencies(examples_run "${target}_run")
set_target_properties("${target}_run" PROPERTIES FOLDER "examples")