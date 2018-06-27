set(sources
  write_matlab.cc
  write_matlab.h
)

add_library(lf.mesh.mesh_utils ${sources})
target_link_libraries(lf.mesh.mesh_utils PUBLIC Eigen3::Eigen Boost::boost lf.mesh)
target_compile_features(lf.mesh.mesh_utils PUBLIC cxx_std_17)
