# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (ellbvp_linfe_demo ellbvp_linfe_demo.cc)

add_executable(lf.examples.ellbvp_linfe.ellbvp_linfe_demo ${ellbvp_linfe_demo})

target_link_libraries(lf.examples.ellbvp_linfe.ellbvp_linfe_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.assemble
  lf.io
  lf.lagrfe
  )

target_compile_features(lf.examples.ellbvp_linfe.ellbvp_linfe_demo PUBLIC cxx_std_17)


