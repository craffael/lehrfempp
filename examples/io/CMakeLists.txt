# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(io_demo io_demo.cc)

add_executable(examples.io.io_demo ${io_demo})

target_link_libraries(examples.io.io_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.io
  )
  
    

  

