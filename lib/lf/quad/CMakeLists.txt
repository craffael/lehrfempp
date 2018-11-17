set(sources
  gauss_quadrature.h
  gauss_quadrature.cc
  make_quad_rule.h
  make_quad_rule.cc
  quad_rule.h
  quad_rule.cc
  quad_rules_tria.cc
)

add_library(quad ${sources})
target_link_libraries(quad PUBLIC Eigen3::Eigen base)

add_subdirectory(test)
