set(sources
  geometry.h
  geometry_interface.h
  point.h
  point.cc
  quad_o1.h
  quad_o1.cc
  segment_o1.h
  segment_o1.cc
  tria_o1.h
  tria_o1.cc
)

add_library(lf.mesh.geometry ${sources})
target_link_libraries(lf.mesh.geometry PUBLIC Eigen3::Eigen base mesh)
target_compile_features(lf.mesh.geometry PUBLIC cxx_std_17)
