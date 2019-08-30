set(comm_sources 
comm_demo.cc 
outside.cc
outside.h
)
add_executable(examples.base.comm.comm_demo ${comm_sources})
target_link_libraries(examples.base.comm.comm_demo PUBLIC Eigen3::Eigen Boost::boost Boost::program_options GTest::main lf.base)
target_compile_features(examples.base.comm.comm_demo PUBLIC cxx_std_17)
