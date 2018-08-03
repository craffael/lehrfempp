# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(refinement_demo refinement_demo.cc)

add_executable(lf.examples.refinement.refinement_demo ${refinement_demo})

target_link_libraries(lf.examples.refinement.refinement_demo
  PUBLIC Eigen3::Eigen Boost::boost GTest::main
  lf.refinement
  lf.mesh.hybrid2dp
  lf.mesh.test_utils
  lf.mesh.utils)

target_compile_features(lf.examples.base.ref_el_output PUBLIC cxx_std_17)

  
    

  

