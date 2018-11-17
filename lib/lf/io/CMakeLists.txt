set(sources
  gmsh_reader.cc
  gmsh_reader.h
  io.h
  vtk_writer.h
  vtk_writer.cc
  write_matplotlib.h
  write_matplotlib.cc
)

lf_add_library(io ${sources})
target_link_libraries(io PUBLIC Eigen3::Eigen base mesh mesh.utils)
if(WIN32) 
  target_compile_options(io PRIVATE "/bigobj")
endif()


add_subdirectory(test)
add_subdirectory(test_utils)
