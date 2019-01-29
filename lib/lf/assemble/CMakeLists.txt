set(sources
  assemble.h
  assembly_types.h
  dofhandler.h
  dofhandler.cc
  coomatrix.h
  coomatrix.cc
  assembler.h
  assembler.cc
  fix_dof.h
  fix_dof.cc
)
lf_add_library(lf.assemble ${sources})
target_link_libraries(lf.assemble PUBLIC
                      Eigen3::Eigen lf.mesh lf.base lf.geometry lf.mesh.utils)

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
