set(sources
  gauss_quadrature.h
  gauss_quadrature.cc
  make_quad_rule.h
  make_quad_rule.cc
  make_quad_rule_nodal.h
  make_quad_rule_nodal.cc
  quad.h
  quad_rule.h
  quad_rule.cc
  quad_rules_tria.cc
)

lf_add_library(lf.quad ${sources})
target_link_libraries(lf.quad PUBLIC Eigen3::Eigen lf.base)
set_target_properties(lf.quad PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
