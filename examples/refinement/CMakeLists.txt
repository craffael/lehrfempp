# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(refinement_demo refinement_demo.cc)

add_executable(examples.refinement.refinement_demo ${refinement_demo})

target_link_libraries(examples.refinement.refinement_demo
        PUBLIC Eigen3::Eigen Boost::boost
        mesh
        mesh.test_utils
        mesh.utils
        mesh.hybrid2d
        refinement
        io
        )


set(point_refinement_demo point_refinement_demo.cc)

add_executable(examples.refinement.point_refinement_demo ${point_refinement_demo})

target_link_libraries(examples.refinement.point_refinement_demo
        PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
        io
        mesh
        mesh.utils
        mesh.hybrid2d
        refinement
        )


set(torus_refinement_demo torus_refinement_demo.cc)

add_executable(examples.refinement.torus_refinement_demo ${torus_refinement_demo})

target_link_libraries(examples.refinement.torus_refinement_demo
        PUBLIC Eigen3::Eigen Boost::filesystem Boost::system Boost::program_options
        io
        mesh
        mesh.hybrid2d
        refinement
        )
