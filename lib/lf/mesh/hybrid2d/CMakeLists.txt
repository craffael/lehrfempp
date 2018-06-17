set(sources
  entity.h
  entity.cc
  hybrid2d.h
  mesh.cc
  mesh.h
  mesh_factory.cc
  mesh_factory.h
)

add_library(lf.mesh.hybrid2d ${sources})
target_link_libraries(lf.mesh.hybrid2d PUBLIC Eigen3::Eigen lf.base lf.mesh)
target_compile_features(lf.mesh.hybrid2d PUBLIC cxx_std_17)

add_subdirectory(test)