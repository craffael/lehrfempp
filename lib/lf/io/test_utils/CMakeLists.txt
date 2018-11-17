
set(sources
  read_mesh.h
  read_mesh.cc
)

add_library(io.test_utils ${sources})
target_link_libraries(io.test_utils PUBLIC Eigen3::Eigen Boost::boost GTest::main io mesh.hybrid2d)
