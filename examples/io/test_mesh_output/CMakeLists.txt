
set(test_mesh_output_sources test_mesh_output_demo.cc)
add_executable(examples.io.test_mesh_output ${test_mesh_output_sources})
target_link_libraries(examples.io.test_mesh_output PUBLIC Eigen3::Eigen
  Boost::boost GTest::main lf.io lf.mesh.hybrid2d lf.mesh.test_utils lf.mesh.utils lf.io)
add_custom_target(examples.io.test_mesh_output_run COMMAND examples.io.test_mesh_output)
add_dependencies(examples_run examples.io.test_mesh_output_run)
