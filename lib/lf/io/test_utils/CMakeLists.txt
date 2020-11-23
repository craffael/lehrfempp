
set(sources
  read_mesh.h
  read_mesh.cc
)

add_library(lf.io.test_utils ${sources})
target_link_libraries(lf.io.test_utils PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main lf.io lf.mesh.hybrid2d)
set_target_properties(lf.io.test_utils PROPERTIES FOLDER "lib")