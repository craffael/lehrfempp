set(comm_sources 
comm_demo.cc 
outside.cc
outside.h
)
add_executable(examples.base.comm.comm_demo ${comm_sources})
target_link_libraries(examples.base.comm.comm_demo PUBLIC Eigen3::Eigen Boost::boost Boost::program_options GTest::gtest_main lf.base)
add_custom_target(examples.base.comm.comm_demo_run COMMAND examples.base.comm.comm_demo)
add_dependencies(examples_run examples.base.comm.comm_demo_run)

