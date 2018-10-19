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

add_library(lf.assemble ${sources})
target_link_libraries(lf.assemble PUBLIC
                      Eigen3::Eigen lf.mesh lf.base lf.geometry lf.mesh.utils)
target_compile_features(lf.assemble PUBLIC cxx_std_17)

add_subdirectory(test)
