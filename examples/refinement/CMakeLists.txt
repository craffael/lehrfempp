# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(refinement_demo refinement_demo.cc)

add_executable(lf.examples.refinement.refinement_demo ${refinement_demo})

target_link_libraries(lf.examples.refinement.refinement_demo
  PUBLIC Eigen3::Eigen Boost::boost
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2dp
  lf.refinement
  lf.io
  )

target_compile_features(lf.examples.refinement.refinement_demo PUBLIC cxx_std_17)


set(point_refinement_demo point_refinement_demo.cc)

add_executable(lf.examples.refinement.point_refinement_demo ${point_refinement_demo})

target_link_libraries(lf.examples.refinement.point_refinement_demo
  PUBLIC Eigen3::Eigen Boost::boost
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2dp
  lf.refinement )

target_compile_features(lf.examples.refinement.point_refinement_demo PUBLIC cxx_std_17)

