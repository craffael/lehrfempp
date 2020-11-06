# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(refinement_demo refinement_demo.cc)

add_executable(examples.refinement.refinement_demo ${refinement_demo})

target_link_libraries(examples.refinement.refinement_demo
        PUBLIC Eigen3::Eigen Boost::boost Boost::program_options
        lf.mesh
        lf.mesh.test_utils
        lf.mesh.utils
        lf.mesh.hybrid2d
        lf.refinement
        lf.io
        )
set_target_properties(examples.refinement.refinement_demo PROPERTIES FOLDER "examples")
add_custom_target(examples.refinement.refinement_demo_run COMMAND examples.refinement.refinement_demo)
add_dependencies(examples_run examples.refinement.refinement_demo_run)
set_target_properties(examples.refinement.refinement_demo_run PROPERTIES FOLDER "examples")


set(point_refinement_demo point_refinement_demo.cc)

add_executable(examples.refinement.point_refinement_demo ${point_refinement_demo})

target_link_libraries(examples.refinement.point_refinement_demo
        PUBLIC Eigen3::Eigen Boost::program_options
        lf.io
        lf.mesh
        lf.mesh.utils
        lf.mesh.hybrid2d
        lf.refinement
        )
set_target_properties(examples.refinement.point_refinement_demo PROPERTIES FOLDER "examples")
add_custom_target(examples.refinement.point_refinement_demo_run COMMAND examples.refinement.point_refinement_demo)
add_dependencies(examples_run examples.refinement.point_refinement_demo_run)
set_target_properties(examples.refinement.point_refinement_demo_run PROPERTIES FOLDER "examples")



set(torus_refinement_demo torus_refinement_demo.cc)

add_executable(examples.refinement.torus_refinement_demo ${torus_refinement_demo})

target_link_libraries(examples.refinement.torus_refinement_demo
        PUBLIC Eigen3::Eigen Boost::program_options
        lf.io
        lf.mesh
        lf.mesh.hybrid2d
        lf.refinement
        )
set_target_properties(examples.refinement.torus_refinement_demo PROPERTIES FOLDER "examples")
add_custom_target(examples.refinement.torus_refinement_demo_run COMMAND examples.refinement.torus_refinement_demo)
add_dependencies(examples_run examples.refinement.torus_refinement_demo_run)
set_target_properties(examples.refinement.torus_refinement_demo_run PROPERTIES FOLDER "examples")



set(second_order_demo second_order_demo.cc)
add_executable(examples.refinement.second_order_demo ${second_order_demo})
target_link_libraries(examples.refinement.second_order_demo
        PUBLIC Eigen3::Eigen
        lf.io
        lf.mesh
        lf.mesh.hybrid2d
        lf.refinement
        )
set_target_properties(examples.refinement.second_order_demo PROPERTIES FOLDER "examples")
add_custom_target(examples.refinement.second_order_demo_run COMMAND examples.refinement.second_order_demo)
add_dependencies(examples_run examples.refinement.second_order_demo_run)
set_target_properties(examples.refinement.second_order_demo_run PROPERTIES FOLDER "examples")
