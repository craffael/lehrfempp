# cmake configuration file for compilation of LehrFEM++ examples dofhandling/assembly

set (dof_demo dof_demo.cc)

add_executable(lf.examples.assemble.dof_demo ${dof_demo})

target_link_libraries(lf.examples.assemble.dof_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.assemble
  )

target_compile_features(lf.examples.assemble.dof_demo PUBLIC cxx_std_17)

set (Dirichlet_Laplacian_demo Dirichlet_Laplacian_demo.cc)

add_executable(lf.examples.assemble.Dirichlet_Laplacian_demo ${Dirichlet_Laplacian_demo})

target_link_libraries(lf.examples.assemble.Dirichlet_Laplacian_demo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.assemble
  lf.io
  lf.fe
  )

target_compile_features(lf.examples.assemble.Dirichlet_Laplacian_demo PUBLIC cxx_std_17)


