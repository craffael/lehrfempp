set(sources
  hybrid2d.h
  mesh.cc
  mesh.h
  mesh_factory.cc
  mesh_factory.h
  point.h
  point.cc
  quad.h
  quad.cc
  segment.h
  segment.cc
  triangle.h
  triangle.cc
)

lf_add_library(mesh.hybrid2d ${sources})
target_link_libraries(mesh.hybrid2d PUBLIC Eigen3::Eigen base mesh mesh.utils)

add_subdirectory(test)
