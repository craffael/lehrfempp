set(sources
  entity.hpp 
  mesh.cpp
  mesh.hpp
)

add_library(mesh ${sources})
target_link_libraries(mesh PUBLIC Eigen3::Eigen)
target_compile_features(mesh PUBLIC cxx_std_17)

