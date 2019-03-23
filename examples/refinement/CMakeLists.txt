# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(refinement_demo refinement_demo.cc)

add_executable(examples.refinement.refinement_demo ${refinement_demo})

target_link_libraries(examples.refinement.refinement_demo
        PUBLIC Eigen3::Eigen Boost::boost
        lf.mesh
        lf.mesh.test_utils
        lf.mesh.utils
        lf.mesh.hybrid2d
        lf.refinement
        lf.io
        )


set(point_refinement_demo point_refinement_demo.cc)

add_executable(examples.refinement.point_refinement_demo ${point_refinement_demo})

target_link_libraries(examples.refinement.point_refinement_demo
        PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
        lf.io
        lf.mesh
        lf.mesh.utils
        lf.mesh.hybrid2d
        lf.refinement
        )


set(torus_refinement_demo torus_refinement_demo.cc)

add_executable(examples.refinement.torus_refinement_demo ${torus_refinement_demo})

target_link_libraries(examples.refinement.torus_refinement_demo
        PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
        lf.io
        lf.mesh
        lf.mesh.hybrid2d
        lf.refinement
        )

set(second_order_demo second_order_demo.cc)
add_executable(examples.refinement.second_order_demo ${second_order_demo})
target_link_libraries(examples.refinement.second_order_demo
        PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
        lf.io
        lf.mesh
        lf.mesh.hybrid2d
        lf.refinement
        )
