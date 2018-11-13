set(sources
  base.h
  comm.cc
  comm.h
  dereference_lambda_random_access_iterator.h
  forward_iterator.h
  forward_range.h
  invalid_type_exception.h
  lf_assert.cc
  lf_assert.h
  lf_exception.h
  random_access_iterator.h
  random_access_range.h
  ref_el.cc
  ref_el.h
  static_vars.h
  static_vars.cc
)

add_library(lf.base ${sources})
target_link_libraries(lf.base PUBLIC Eigen3::Eigen Boost::boost Boost::program_options)
target_compile_features(lf.base PUBLIC cxx_std_17)
target_include_directories(lf.base PUBLIC
  "$<BUILD_INTERFACE:${LOCAL_INCLUDE_DIRECTORY}>"
  "$<BUILD_INTERFACE:${LOCAL_GENERATED_INCLUDE_DIRECTORY}>"
)

add_subdirectory(test)
