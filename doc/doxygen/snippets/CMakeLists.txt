# Here we want to compile all snippets to make sure that they stay up-to-date, even when the code is refactored.

set(sources 
  ref_el.cc
)

add_library(docsnippets ${sources})
get_target_property(OUT base INCLUDE_DIRECTORIES)
message("out = ${OUT}")
target_link_libraries(docsnippets PUBLIC Eigen3::Eigen base)
target_compile_features(docsnippets PUBLIC cxx_std_17)

message("hello ${LOCAL_INCLUDE_DIRECTORY}")