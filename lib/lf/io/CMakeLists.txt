set(sources
  gmsh_reader.cc
  gmsh_reader.h
  io.h
  vtk_writer.h
  vtk_writer.cc
  write_matplotlib.h
  write_matplotlib.cc
)

lf_add_library(lf.io ${sources})
target_link_libraries(lf.io PUBLIC Eigen3::Eigen lf.base lf.mesh lf.mesh.utils)
if(WIN32 AND NOT MINGW) 
  target_compile_options(lf.io PRIVATE "/bigobj")
endif()

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
add_subdirectory(test_utils)
