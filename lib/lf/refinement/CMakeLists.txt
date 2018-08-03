set(sources
  refinement.h refinement.cc
  refinement_hierarchy.h refinement_hierarchy.cc
  refutils.h refutils.cc
  )

add_library(lf.refinement ${sources})
target_link_libraries(lf.refinement PUBLIC Eigen3::Eigen lf.base lf.geometry)
target_compile_features(lf.refinement PUBLIC cxx_std_17)

add_subdirectory(test)
