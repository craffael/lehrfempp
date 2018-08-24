set(sources
  gauss_legendre.h
  gauss_legendre.cc
  make_quad_rule.h
  make_quad_rule.cc
  quad_rule.h
  quad_rule.cc
)

add_library(lf.quad ${sources})
target_link_libraries(lf.quad PUBLIC Eigen3::Eigen lf.base)
target_compile_features(lf.mesh PUBLIC cxx_std_17)

add_subdirectory(test)
