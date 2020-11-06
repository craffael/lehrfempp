set(sources
  gmsh_file_v2.cc
  gmsh_file_v2.h
  gmsh_file_v4.cc
  gmsh_file_v4.h
  gmsh_file_v4_detail.h
  gmsh_file_v4_text.cc
  gmsh_file_v4_binary.cc
  gmsh_reader.cc
  gmsh_reader.h
  io.h
  vtk_writer.h
  vtk_writer.cc
  write_matlab.h
  write_matlab.cc
  write_matplotlib.h
  write_matplotlib.cc
  write_tikz.h
  write_tikz.cc
)

lf_add_library(lf.io ${sources})
target_link_libraries(lf.io PUBLIC Eigen3::Eigen lf.base lf.mesh lf.mesh.utils)
set_target_properties(lf.io PROPERTIES FOLDER "lib")
if(WIN32 AND NOT MINGW) 
  target_compile_options(lf.io PRIVATE "/bigobj")
endif()

if(LF_ENABLE_TESTING)
  add_subdirectory(test)
endif()
add_subdirectory(test_utils)
