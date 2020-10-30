# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (dof_demo dof_demo.cc)

add_executable(examples.assemble.dof_demo ${dof_demo})
target_link_libraries(examples.assemble.dof_demo
  PUBLIC Eigen3::Eigen Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.assemble
  )
add_custom_target(examples.assemble.dof_demo_run COMMAND examples.assemble.dof_demo)
add_dependencies(examples_run examples.assemble.dof_demo_run)



set (Dirichlet_Laplacian_demo Dirichlet_Laplacian_demo.cc)
add_executable(examples.assemble.Dirichlet_Laplacian_demo ${Dirichlet_Laplacian_demo})
target_link_libraries(examples.assemble.Dirichlet_Laplacian_demo
  PUBLIC Eigen3::Eigen Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.assemble
  lf.io
  lf.uscalfe
  )
add_custom_target(examples.assemble.Dirichlet_Laplacian_demo_run COMMAND examples.assemble.Dirichlet_Laplacian_demo)
add_dependencies(examples_run examples.assemble.Dirichlet_Laplacian_demo_run)


