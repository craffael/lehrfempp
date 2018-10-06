set(sources
  all_codim_mesh_data_set.h
  codim_mesh_data_set.h
  lambda_mesh_data_set.h
  mesh_data_set.h
  print_info.cc
  print_info.h
  utils.h
  write_matlab.cc
  write_matlab.h
  write_matplotlib.cc
  write_matplotlib.h
  write_tikz.cc
  write_tikz.h
)

add_library(lf.mesh.utils ${sources})
target_link_libraries(lf.mesh.utils PUBLIC Eigen3::Eigen Boost::boost lf.mesh)
target_compile_features(lf.mesh.utils PUBLIC cxx_std_17)
