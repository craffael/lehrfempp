# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (ellbvp_msh_demo ellbvp_msh_demo.cc)

add_executable(examples.ellbvp_msh.ellbvp_msh_demo ${ellbvp_msh_demo})

target_link_libraries(examples.ellbvp_msh.ellbvp_msh_demo
  PUBLIC Eigen3::Eigen Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.assemble
  lf.io
  lf.uscalfe
  )

add_custom_target(examples.ellbvp_msh.ellbvp_msh_demo_run COMMAND examples.ellbvp_msh.ellbvp_msh_demo)
add_dependencies(examples_run examples.ellbvp_msh.ellbvp_msh_demo_run)


