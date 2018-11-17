# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(io_demo io_demo.cc)

add_executable(examples.io.io_demo ${io_demo})

target_link_libraries(examples.io.io_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  mesh
  mesh.test_utils
  mesh.utils
  mesh.hybrid2d
  refinement
  io
  )
  
    

  

