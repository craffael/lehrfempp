set(sources
  refinement.h refinement.cc
  hybrid2d_refinement_pattern.h hybrid2d_refinement_pattern.cc
  mesh_hierarchy.h mesh_hierarchy.cc
  refutils.h refutils.cc
  )

add_library(lf.refinement ${sources})
target_link_libraries(lf.refinement PUBLIC Eigen3::Eigen lf.base lf.geometry)
target_compile_features(lf.refinement PUBLIC cxx_std_17)

add_subdirectory(test)
