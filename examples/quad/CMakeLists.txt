# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(quad_demo quad_demo.cc)

add_executable(lf.examples.quad.quad_demo ${quad_demo})

target_link_libraries(lf.examples.quad.quad_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2dp
  lf.refinement
  lf.io
  lf.quad
  )

target_compile_features(lf.examples.quad.quad_demo PUBLIC cxx_std_17)
  
    

  

