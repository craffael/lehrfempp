# Set sources
set(SOURCES 
  strangsplitting.h
  strangsplitting.cc)

# Define executable
add_executable(splitting strangsplitting_main.cc ${SOURCES})
add_executable(testIsland modelproblem_island_main.cc  ${SOURCES})
add_executable(testCircle modelproblem_circle_main.cc  ${SOURCES})
add_executable(test3Circles modelproblem_threecircles_main.cc  ${SOURCES})
add_executable(convergenceH convergence_meshwidth_main.cc  ${SOURCES})
add_executable(convergenceTau convergence_tau_main.cc  ${SOURCES})
add_executable(convergenceLinked convergence_linked_main.cc  ${SOURCES})

# link executable with required modules

target_link_libraries(splitting
  PUBLIC Boost::program_options  Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

target_link_libraries(testIsland
  PUBLIC Boost::program_options  Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

target_link_libraries(testCircle
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

target_link_libraries(test3Circles
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

target_link_libraries(convergenceH
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

target_link_libraries(convergenceTau
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

target_link_libraries(convergenceLinked
  PUBLIC Boost::program_options Eigen3::Eigen
  lf.assemble
  lf.base
  lf.geometry
  lf.io
  lf.mesh
  lf.mesh.hybrid2d
  lf.mesh.test_utils
  lf.mesh.utils
  lf.refinement
  lf.uscalfe
  )

add_custom_target(splitting_run COMMAND splitting)
add_custom_target(testIsland_run COMMAND testIsland)
add_custom_target(testCircle_run COMMAND testCircle)
add_custom_target(test3Circles_run COMMAND test3Circles)
add_custom_target(convergenceH_run COMMAND convergenceH)
add_custom_target(convergenceTau_run COMMAND convergenceTau)
add_custom_target(convergenceLinked_run COMMAND convergenceLinked)
