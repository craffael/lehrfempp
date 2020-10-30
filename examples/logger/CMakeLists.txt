

add_executable(examples.logger.mesh_hierarchy_demo mesh_hierarchy_demo.cc)
target_link_libraries(examples.logger.mesh_hierarchy_demo
  PUBLIC Eigen3::Eigen Boost::system Boost::program_options
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.io
  )
add_custom_target(examples.logger.mesh_hierarchy_demo_run COMMAND examples.assemble.dof_demo)
add_dependencies(examples_run examples.logger.mesh_hierarchy_demo_run)





