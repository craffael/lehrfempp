
set(sources
  test_entity_indexing.cc
  test_entity_indexing.h
)

add_library(lf.mesh.test_utils ${sources})
target_link_libraries(lf.mesh.test_utils PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh)
target_compile_features(lf.mesh.test_utils PUBLIC cxx_std_17)
