set(src
  geom.h
  brep_lagr_segment.h
  brep_lagr_segment.cc
  brep_lagr_tria.h
  brep_curve.h
  brep_curve.cc
  brep_surface_segment.h
  brep_surface_segment.cc
  brep_tria_transfinite.h
  brep_tria_transfinite.cc
)

add_library(lf.brep.geom ${src})
target_link_libraries(lf.brep.geom PUBLIC lf.base lf.brep.occt)
set_target_properties(lf.brep.geom PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()