set(comm_sources comm_demo.cc outside.cc)
add_executable(lf.examples.base.comm.comm_demo ${comm_sources})
target_link_libraries(lf.examples.base.comm.comm_demo PUBLIC Eigen3::Eigen Boost::boost Boost::program_options GTest::main lf.mesh.test_utils lf.mesh.utils lf.geometry lf.refinement)
target_compile_features(lf.examples.base.comm.comm_demo PUBLIC cxx_std_17)
