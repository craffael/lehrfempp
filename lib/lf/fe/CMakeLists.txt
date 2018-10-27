set(sources
  lagr_fe.h
  lagr_fe.cc
)

add_library(lf.fe ${sources})
target_link_libraries(lf.fe PUBLIC
  Eigen3::Eigen lf.mesh lf.base lf.geometry
  lf.mesh.utils lf.assemble)
target_compile_features(lf.fe PUBLIC cxx_std_17)

# add_subdirectory(test)
