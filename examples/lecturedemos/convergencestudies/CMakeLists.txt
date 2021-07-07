add_executable(examples.lecturedemos.convergencestudies.linquadlagrfe linquadlagrfe.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.linquadlagrfe
    PUBLIC Eigen3::Eigen
    Boost::program_options
    lf.base
    lf.uscalfe
    lf.mesh
    lf.mesh.utils
    lf.mesh.hybrid2d
    lf.assemble
    lf.quad
    lf.io
)
set_target_properties(examples.lecturedemos.convergencestudies.linquadlagrfe PROPERTIES FOLDER "examples")
add_custom_target(examples.lecturedemos.convergencestudies.linquadlagrfe_run
    COMMAND examples.lecturedemos.convergencestudies.linquadlagrfe -o ${CMAKE_CURRENT_BINARY_DIR}/results_linquadlagrfe.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/linquadlagrfe_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results_linquadlagrfe.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linquadlagrfe_plot.py
)
add_dependencies(examples_run examples.lecturedemos.convergencestudies.linquadlagrfe_run)
set_target_properties(examples_run examples.lecturedemos.convergencestudies.linquadlagrfe_run PROPERTIES FOLDER "examples")

add_executable(examples.lecturedemos.convergencestudies.linquadfelshaped linquadfelshaped.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.linquadfelshaped
    PUBLIC Eigen3::Eigen
    Boost::program_options
    lf.base
    lf.uscalfe
    lf.mesh
    lf.mesh.utils
    lf.mesh.hybrid2d
    lf.assemble
    lf.quad
    lf.io
)
set_target_properties(examples.lecturedemos.convergencestudies.linquadfelshaped PROPERTIES FOLDER "examples")
add_custom_target(examples.lecturedemos.convergencestudies.linquadfelshaped_run
    COMMAND examples.lecturedemos.convergencestudies.linquadfelshaped -o ${CMAKE_CURRENT_BINARY_DIR}/results_linquadfelshaped.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/linquadfelshaped_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results_linquadfelshaped.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linquadfelshaped_plot.py
)
add_dependencies(examples_run examples.lecturedemos.convergencestudies.linquadfelshaped_run)
set_target_properties(examples_run examples.lecturedemos.convergencestudies.linquadfelshaped_run PROPERTIES FOLDER "examples")

add_executable(examples.lecturedemos.convergencestudies.prefinement prefinement.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.prefinement
    PUBLIC Eigen3::Eigen
    Boost::program_options
    lf.base
    lf.uscalfe
    lf.mesh
    lf.mesh.utils
    lf.mesh.hybrid2d
    lf.assemble
    lf.quad
    lf.io
    lf.refinement
)
set_target_properties(examples.lecturedemos.convergencestudies.prefinement PROPERTIES FOLDER "examples")
add_custom_target(examples.lecturedemos.convergencestudies.prefinement_run
    COMMAND examples.lecturedemos.convergencestudies.prefinement -o ${CMAKE_CURRENT_BINARY_DIR}/results_prefinement.csv -p 19 && python3 ${CMAKE_CURRENT_SOURCE_DIR}/prefinement_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results_prefinement.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/prefinement_plot.py
)
add_dependencies(examples_run examples.lecturedemos.convergencestudies.prefinement_run)
set_target_properties(examples_run examples.lecturedemos.convergencestudies.prefinement_run PROPERTIES FOLDER "examples")
