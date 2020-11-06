set(sources
  entity.h
  entity.cc
  mesh.h
  mesh.cc
  mesh_factory.h
  mesh_interface.h
  mesh_interface.cc
)

lf_add_library(lf.mesh ${sources})
target_link_libraries(lf.mesh PUBLIC Eigen3::Eigen lf.base lf.geometry lf.mesh.utils)
set_target_properties(lf.mesh PROPERTIES FOLDER "lib")

add_subdirectory(hybrid2d)
add_subdirectory(utils)
add_subdirectory(test_utils)
