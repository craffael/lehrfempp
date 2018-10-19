# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (dof_demo dof_demo.cc)

add_executable(lf.examples.assemble.dof_demo ${dof_demo})

target_link_libraries(lf.examples.assemble.dof_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.hybrid2dp
  lf.refinement
  lf.assemble
  )

target_compile_features(lf.examples.io.io_demo PUBLIC cxx_std_17)


