set(sources
  print_info.cc
  print_info.h
  utils.h
  write_matlab.cc
  write_matlab.h
)

add_library(lf.mesh.utils ${sources})
target_link_libraries(lf.mesh.utils PUBLIC Eigen3::Eigen Boost::boost lf.mesh)
target_compile_features(lf.mesh.utils PUBLIC cxx_std_17)
