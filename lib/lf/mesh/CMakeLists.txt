set(sources
  entity.h
  mesh.cc
  mesh.h
  mesh_builder.h
  mesh_interface.h
)

add_library(lf.mesh ${sources})
target_link_libraries(lf.mesh PUBLIC Eigen3::Eigen lf.base)
target_compile_features(lf.mesh PUBLIC cxx_std_17)

add_subdirectory(geometry)
add_subdirectory(hybrid2d)
