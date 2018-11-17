# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(quad_demo quad_demo.cc)

add_executable(examples.quad.quad_demo ${quad_demo})

target_link_libraries(examples.quad.quad_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  mesh
  mesh.test_utils
  mesh.utils
  mesh.hybrid2d
  refinement
  io
  quad
  )
  
    

  

