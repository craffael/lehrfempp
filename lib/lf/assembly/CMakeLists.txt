set(sources
  dofhandler.h
  dofhandler.cc
  local_static_dofs.h
  local_static_dofs.cc
)

add_library(lf.assembly ${sources})
target_link_libraries(lf.assembly PUBLIC
                      Eigen3::Eigen lf.mesh lf.base lf.geometry lf.mesh.utils)
target_compile_features(lf.assembly PUBLIC cxx_std_17)

add_subdirectory(test)
