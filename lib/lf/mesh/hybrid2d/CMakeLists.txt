set(sources
  entity.h
  entity.cc
  geometry.h
  geometry.cc
  mesh.cc
  mesh.h
)

add_library(lf.mesh.hybrid2d ${sources})
target_link_libraries(lf.mesh.hybrid2d PUBLIC Eigen3::Eigen base mesh)
target_compile_features(lf.mesh.hybrid2d PUBLIC cxx_std_17)
