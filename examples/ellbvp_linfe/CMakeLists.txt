# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (ellbvp_linfe_demo ellbvp_linfe_demo.cc)

add_executable(examples.ellbvp_linfe.ellbvp_linfe_demo ${ellbvp_linfe_demo})

target_link_libraries(examples.ellbvp_linfe.ellbvp_linfe_demo
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
set_target_properties(examples.ellbvp_linfe.ellbvp_linfe_demo PROPERTIES FOLDER "examples")

add_custom_target(examples.ellbvp_linfe.ellbvp_linfe_demo_run COMMAND examples.ellbvp_linfe.ellbvp_linfe_demo)
add_dependencies(examples_run examples.ellbvp_linfe.ellbvp_linfe_demo_run)
set_target_properties(examples.ellbvp_linfe.ellbvp_linfe_demo_run PROPERTIES FOLDER "examples")


