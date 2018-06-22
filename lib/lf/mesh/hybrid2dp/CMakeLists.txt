set(sources
  hybrid2dp.h
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

add_library(lf.mesh.hybrid2dp ${sources})
target_link_libraries(lf.mesh.hybrid2dp PUBLIC Eigen3::Eigen lf.base lf.mesh)
target_compile_features(lf.mesh.hybrid2dp PUBLIC cxx_std_17)

#add_subdirectory(test)
