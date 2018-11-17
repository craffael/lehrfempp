set(sources
  refinement.h refinement.cc
  hybrid2d_refinement_pattern.h hybrid2d_refinement_pattern.cc
  mesh_hierarchy.h mesh_hierarchy.cc
  refutils.h refutils.cc
  )

add_library(refinement ${sources})
target_link_libraries(refinement PUBLIC Eigen3::Eigen base geometry)

add_subdirectory(test)
