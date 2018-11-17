set(sources
  entity.h
  entity.cc
  mesh.h
  mesh.cc
  mesh_factory.h
  mesh_interface.h
  mesh_interface.cc
)

lf_add_library(mesh ${sources})
target_link_libraries(mesh PUBLIC Eigen3::Eigen base geometry mesh.utils)

add_subdirectory(hybrid2d)
add_subdirectory(utils)
add_subdirectory(test_utils)
