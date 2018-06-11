set(sources
  entity.h
  mesh.h
  mesh.cc
  mesh_builder.h
  mesh_interface.h
)

add_library(lf.mesh ${sources})
target_link_libraries(lf.mesh PUBLIC Eigen3::Eigen lf.base lf.geometry)
target_compile_features(lf.mesh PUBLIC cxx_std_17)

add_subdirectory(hybrid2d)
add_subdirectory(test_utils)
