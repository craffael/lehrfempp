set(sources
  gmsh_reader.cc
  gmsh_reader.h
  io.h
)

add_library(lf.io ${sources})
target_link_libraries(lf.io PUBLIC Eigen3::Eigen lf.base lf.mesh)
target_compile_features(lf.io PUBLIC cxx_std_17)
if(WIN32) 
  target_compile_options(lf.io PRIVATE "/bigobj")
endif()


add_subdirectory(test)