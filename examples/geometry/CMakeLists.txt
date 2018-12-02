#cmake configuration file for compilation of LehrFEM++ examples: Geometry

set(geometry_sources geometry_output.cc)

add_executable(examples.geometry.geometry_output ${geometry_sources})

# Target link libraries? See similar file for base example
target_link_libraries(examples.geometry.geometry_output PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.hybrid2d lf.mesh.test_utils lf.mesh.utils)

