set(sources
  assemble.h
  assembly_types.h
  dofhandler.h
  dofhandler.cc
  coomatrix.h
  coomatrix.cc
  assembler.h
  assembler.cc
)
lf_add_library(assemble ${sources})
target_link_libraries(assemble PUBLIC
                      Eigen3::Eigen mesh base geometry mesh.utils)

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
