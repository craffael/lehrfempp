set(sources 
  base.h
  forward_iterator.h
  lf_assert.h
  forward_range.h
  ref_el.cc
  ref_el.h
)



add_library(base ${sources})
target_link_libraries(base PUBLIC Eigen3::Eigen Boost::boost)
target_compile_features(base PUBLIC cxx_std_17)
target_include_directories(base PUBLIC
  "$<BUILD_INTERFACE:${LOCAL_INCLUDE_DIRECTORY}>"
  "$<BUILD_INTERFACE:${LOCAL_GENERATED_INCLUDE_DIRECTORY}>"
)




add_subdirectory(test)