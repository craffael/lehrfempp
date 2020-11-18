
set(interface_src
  brep_model.h
  brep_geometry.h
  interface.h
)

add_library(lf.brep.interface ${interface_src})
target_link_libraries(lf.brep.interface PUBLIC lf.base)
set_target_properties(lf.brep.interface PROPERTIES FOLDER "lib")
