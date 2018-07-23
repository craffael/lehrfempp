
set(sources
  read_mesh.h
  read_mesh.cc
)

add_library(lf.io.test_utils ${sources})
target_link_libraries(lf.io.test_utils PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.io)
target_compile_features(lf.io.test_utils PUBLIC cxx_std_17)
