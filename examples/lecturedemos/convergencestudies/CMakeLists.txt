

add_executable(examples.lecturedemos.convergencestudies.twopointbvp twopointbvp.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.twopointbvp
    PUBLIC Eigen3::Eigen
    Boost::program_options
    lf.quad
)
add_custom_target(examples.lecturedemos.convergencestudies.twopointbvp_run
    COMMAND examples.lecturedemos.convergencestudies.twopointbvp -o ${CMAKE_CURRENT_BINARY_DIR}/results.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/twopointbvp_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/twopointbvp_plot.py
)
add_dependencies(examples_run examples.lecturedemos.convergencestudies.twopointbvp_run)


add_executable(examples.lecturedemos.convergencestudies.linquadlagrfe linquadlagrfe.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.linquadlagrfe
    PUBLIC Eigen3::Eigen
    Boost::filesystem
    Boost::system
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
add_custom_target(examples.lecturedemos.convergencestudies.linquadlagrfe_run
    COMMAND examples.lecturedemos.convergencestudies.linquadlagrfe -o ${CMAKE_CURRENT_BINARY_DIR}/results.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/linquadlagrfe_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linquadlagrfe_plot.py
)
add_dependencies(examples_run examples.lecturedemos.convergencestudies.linquadlagrfe_run)


add_executable(examples.lecturedemos.convergencestudies.linquadfelshaped linquadfelshaped.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.linquadfelshaped
    PUBLIC Eigen3::Eigen
    Boost::filesystem
    Boost::system
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
add_custom_target(examples.lecturedemos.convergencestudies.linquadfelshaped_run
    COMMAND examples.lecturedemos.convergencestudies.linquadfelshaped -o ${CMAKE_CURRENT_BINARY_DIR}/results.csv && python3 ${CMAKE_CURRENT_SOURCE_DIR}/linquadfelshaped_plot.py ${CMAKE_CURRENT_BINARY_DIR}/results.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/linquadfelshaped_plot.py
)
add_dependencies(examples_run examples.lecturedemos.convergencestudies.linquadfelshaped_run)
