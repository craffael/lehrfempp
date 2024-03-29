set(sources
  geometry.h
  geometry.cc
  geometry_interface.h
  geometry_interface.cc
  point.h
  point.cc
  quad_o1.h
  quad_o1.cc
  quad_o2.h
  quad_o2.cc
  refinement_pattern.h
  refinement_pattern.cc
  segment_o1.h
  segment_o1.cc
  segment_o2.h
  segment_o2.cc
  tria_o1.h
  tria_o1.cc
  tria_o2.h
  tria_o2.cc
  print_info.h
  print_info.cc
)

lf_add_library(lf.geometry ${sources})
target_link_libraries(lf.geometry PUBLIC Eigen3::Eigen lf.base)
set_target_properties(lf.geometry PROPERTIES FOLDER "lib")

if(MSVC)
  target_compile_options(lf.geometry PRIVATE /bigobj)
endif()

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()

add_subdirectory(test_utils)
