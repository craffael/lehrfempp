
set(sources
  check_entity_indexing.cc
  check_entity_indexing.h
  check_geometry_orientation.cc
  check_geometry_orientation.h
  check_local_topology.h
  check_local_topology.cc
  check_mesh_completeness.h
  check_mesh_completeness.cc
  test_meshes.h
  test_meshes.cc
)

lf_add_library(lf.mesh.test_utils ${sources})
target_link_libraries(lf.mesh.test_utils PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh lf.mesh.hybrid2d)
