set(occt_src
  occt.h
  occt_brep_model.h
  occt_brep_model.cc
  occt_curve_geometry.h
  occt_curve_geometry.cc
  occt_utils.h
  occt_utils.cc
)

add_library(lf.brep.occt ${occt_src})
target_link_libraries(lf.brep.occt PUBLIC lf.base  occt::TKTopAlgo  occt::TKGeomAlgo occt::TKBRep occt::TKGeomBase occt::TKG3d  occt::TKG2d  occt::TKMath occt::TKernel)
set_target_properties(lf.brep.occt PROPERTIES FOLDER "lib")

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()