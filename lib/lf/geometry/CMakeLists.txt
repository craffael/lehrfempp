set(sources
  geometry.h
  geometry.cc
  geometry_interface.h
  geometry_interface.cc
  point.h
  point.cc
  quad_o1.h
  quad_o1.cc
  refinement_pattern.h
  refinement_pattern.cc
  segment_o1.h
  segment_o1.cc
  tria_o1.h
  tria_o1.cc
  print_info.h
  print_info.cc
)

add_library(lf.geometry ${sources})
target_link_libraries(lf.geometry PUBLIC Eigen3::Eigen lf.base)
target_compile_features(lf.geometry PUBLIC cxx_std_17)

add_subdirectory(test)
