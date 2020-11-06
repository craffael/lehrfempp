# Set sources
set(SOURCES 
  strangsplitting.h
  strangsplitting.cc)

# Define executable
add_executable(projects.fisher_kpp.splitting strangsplitting_main.cc ${SOURCES})
add_executable(projects.fisher_kpp.testIsland modelproblem_island_main.cc  ${SOURCES})
add_executable(projects.fisher_kpp.testCircle modelproblem_circle_main.cc  ${SOURCES})
add_executable(projects.fisher_kpp.test3Circles modelproblem_threecircles_main.cc  ${SOURCES})
add_executable(projects.fisher_kpp.convergenceH convergence_meshwidth_main.cc  ${SOURCES})
add_executable(projects.fisher_kpp.convergenceTau convergence_tau_main.cc  ${SOURCES})
add_executable(projects.fisher_kpp.convergenceLinked convergence_linked_main.cc  ${SOURCES})

# link executable with required modules

target_link_libraries(projects.fisher_kpp.splitting
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

target_link_libraries(projects.fisher_kpp.testIsland
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

target_link_libraries(projects.fisher_kpp.testCircle
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

target_link_libraries(projects.fisher_kpp.test3Circles
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

target_link_libraries(projects.fisher_kpp.convergenceH
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

target_link_libraries(projects.fisher_kpp.convergenceTau
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

target_link_libraries(projects.fisher_kpp.convergenceLinked
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

set_target_properties(projects.fisher_kpp.splitting         PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.testIsland        PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.testCircle        PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.test3Circles      PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.convergenceH      PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.convergenceTau    PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.convergenceLinked PROPERTIES FOLDER "projects/fisher_kpp")



add_custom_target(projects.fisher_kpp.splitting_run         COMMAND splitting)
add_custom_target(projects.fisher_kpp.testIsland_run        COMMAND testIsland)
add_custom_target(projects.fisher_kpp.testCircle_run        COMMAND testCircle)
add_custom_target(projects.fisher_kpp.test3Circles_run      COMMAND test3Circles)
add_custom_target(projects.fisher_kpp.convergenceH_run      COMMAND convergenceH)
add_custom_target(projects.fisher_kpp.convergenceTau_run    COMMAND convergenceTau)
add_custom_target(projects.fisher_kpp.convergenceLinked_run COMMAND convergenceLinked)

set_target_properties(projects.fisher_kpp.splitting_run         PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.testIsland_run        PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.testCircle_run        PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.test3Circles_run      PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.convergenceH_run      PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.convergenceTau_run    PROPERTIES FOLDER "projects/fisher_kpp")
set_target_properties(projects.fisher_kpp.convergenceLinked_run PROPERTIES FOLDER "projects/fisher_kpp")
