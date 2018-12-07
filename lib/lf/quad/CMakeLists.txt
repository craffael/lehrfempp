set(sources
  gauss_quadrature.h
  gauss_quadrature.cc
  make_quad_rule.h
  make_quad_rule.cc
  quad_rule.h
  quad_rule.cc
  quad_rules_tria.cc
)

add_library(lf.quad ${sources})
target_link_libraries(lf.quad PUBLIC Eigen3::Eigen lf.base)

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
