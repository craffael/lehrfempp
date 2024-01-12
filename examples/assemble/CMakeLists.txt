# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (dof_demo dof_demo.cc)

add_executable(examples.assemble.dof_demo ${dof_demo})
target_link_libraries(examples.assemble.dof_demo
  PUBLIC Eigen3::Eigen Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.assemble
  )
set_target_properties(examples.assemble.dof_demo PROPERTIES FOLDER "examples")
add_custom_target(examples.assemble.dof_demo_run COMMAND examples.assemble.dof_demo)
add_dependencies(examples_run examples.assemble.dof_demo_run)
set_target_properties(examples.assemble.dof_demo_run PROPERTIES FOLDER "examples")

set (dofh_probe dofh_probe.cc)
add_executable(examples.assemble.dofh_probe ${dofh_probe})
target_link_libraries(examples.assemble.dofh_probe
  PUBLIC Eigen3::Eigen Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.assemble
  )
set_target_properties(examples.assemble.dofh_probe PROPERTIES FOLDER "examples")
add_custom_target(examples.assemble.dofh_probe_run COMMAND examples.assemble.dofh_probe)
add_dependencies(examples_run examples.assemble.dofh_probe_run)
set_target_properties(examples.assemble.dofh_probe_run PROPERTIES FOLDER "examples")

set (Dirichlet_Laplacian_demo Dirichlet_Laplacian_demo.cc)
add_executable(examples.assemble.Dirichlet_Laplacian_demo ${Dirichlet_Laplacian_demo})
target_link_libraries(examples.assemble.Dirichlet_Laplacian_demo
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
set_target_properties(examples.assemble.Dirichlet_Laplacian_demo PROPERTIES FOLDER "examples")
add_custom_target(examples.assemble.Dirichlet_Laplacian_demo_run COMMAND examples.assemble.Dirichlet_Laplacian_demo)
add_dependencies(examples_run examples.assemble.Dirichlet_Laplacian_demo_run)
set_target_properties(examples.assemble.Dirichlet_Laplacian_demo_run PROPERTIES FOLDER "examples")

