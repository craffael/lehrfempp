set(sources
  gmsh_reader.cc
  gmsh_reader.h
)

add_library(lf.io ${sources})
target_link_libraries(lf.io PUBLIC Eigen3::Eigen lf.base lf.mesh)
target_compile_features(lf.io PUBLIC cxx_std_17)

add_subdirectory(test)