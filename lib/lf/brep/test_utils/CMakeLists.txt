set(src
  check_brep_geometry.h
  check_brep_geometry.cc
  curve_circle.h
  curve_circle.cc
  
)

add_library(lf.brep.test_utils ${src})
target_link_libraries(lf.brep.test_utils PUBLIC lf.base lf.brep.interface)
set_target_properties(lf.brep.test_utils PROPERTIES FOLDER "lib")


if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()