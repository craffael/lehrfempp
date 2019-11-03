#cmake configuration file for compilation of LehrFEM++ examples: Geometry

set(geometry_output_sources geometry_output.cc)

add_executable(examples.geometry.geometry_output ${geometry_output_sources})

# Target link libraries? See similar file for base example
target_link_libraries(examples.geometry.geometry_output PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main lf.mesh.hybrid2d lf.mesh.test_utils lf.mesh.utils)
add_custom_target(examples.geometry.geometry_output_run COMMAND examples.geometry.geometry_output)
add_dependencies(examples_run examples.geometry.geometry_output_run)
