

add_executable(projects.dpg.convergence_studies.convergence_study1
	       convergence_study1.cc)
add_executable(projects.dpg.convergence_studies.convergence_study2
	       convergence_study2.cc)
add_executable(projects.dpg.convergence_studies.convergence_study3
	       convergence_study3.cc)
add_executable(projects.dpg.convergence_studies.convergence_study4
	       convergence_study4.cc)
add_executable(projects.dpg.convergence_studies.convergence_study5
	       convergence_study5.cc)


target_link_libraries(projects.dpg.convergence_studies.convergence_study1
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.quad
  lf.assemble
  lf.io
  lf.uscalfe
  projects.dpg
  )

target_link_libraries(projects.dpg.convergence_studies.convergence_study2
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.quad
  lf.assemble
  lf.io
  lf.uscalfe
  projects.dpg
  )

target_link_libraries(projects.dpg.convergence_studies.convergence_study3
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.quad
  lf.assemble
  lf.io
  lf.uscalfe
  projects.dpg
  )

target_link_libraries(projects.dpg.convergence_studies.convergence_study4
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.quad
  lf.assemble
  lf.io
  lf.uscalfe
  projects.dpg
  )

target_link_libraries(projects.dpg.convergence_studies.convergence_study5
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
  lf.mesh
  lf.mesh.utils
  lf.mesh.test_utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.quad
  lf.assemble
  lf.io
  lf.uscalfe
  projects.dpg
  )


target_compile_features(projects.dpg.convergence_studies.convergence_study1 PUBLIC cxx_std_17)
target_compile_features(projects.dpg.convergence_studies.convergence_study2 PUBLIC cxx_std_17)
target_compile_features(projects.dpg.convergence_studies.convergence_study3 PUBLIC cxx_std_17)
target_compile_features(projects.dpg.convergence_studies.convergence_study4 PUBLIC cxx_std_17)
target_compile_features(projects.dpg.convergence_studies.convergence_study5 PUBLIC cxx_std_17)

if(WIN32)
    target_compile_options(projects.dpg.convergence_studies.convergence_study1 PRIVATE /bigobj)
    target_compile_options(projects.dpg.convergence_studies.convergence_study2 PRIVATE /bigobj)
    target_compile_options(projects.dpg.convergence_studies.convergence_study3 PRIVATE /bigobj)
    target_compile_options(projects.dpg.convergence_studies.convergence_study4 PRIVATE /bigobj)
    target_compile_options(projects.dpg.convergence_studies.convergence_study5 PRIVATE /bigobj)
  endif()