set(sources 
  base.h
  ref_el.cc
  ref_el.h
  lf_assert.h
)

add_library(base ${sources})
target_link_libraries(base PUBLIC Eigen3::Eigen Boost::boost)
target_compile_features(base PUBLIC cxx_std_17)
message("local include directory: ${LOCAL_INCLUDE_DIRECTORY}")
target_include_directories(base PUBLIC
  "$<BUILD_INTERFACE:${LOCAL_INCLUDE_DIRECTORY}>"
  "$<BUILD_INTERFACE:${LOCAL_GENERATED_INCLUDE_DIRECTORY}>"
)

