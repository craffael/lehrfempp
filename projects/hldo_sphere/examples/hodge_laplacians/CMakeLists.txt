project(projects.hldo_sphere.examples.hodge_laplacians)

add_executable(${PROJECT_NAME}.main hodge_laplacians.cc)
target_link_libraries(${PROJECT_NAME}.main PUBLIC Eigen3::Eigen projects.hldo_sphere.discretization lf.base lf.io lf.quad lf.mesh lf.mesh.hybrid2d lf.refinement projects.hldo_sphere.mesh projects.hldo_sphere.post_processing projects.hldo_sphere.discretization)
set_target_properties(${PROJECT_NAME}.main PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(${PROJECT_NAME}.main PUBLIC cxx_std_17)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output_regular.dat
		   DEPENDS ${PROJECT_NAME}.main
		   COMMAND $<TARGET_FILE:${PROJECT_NAME}.main> regular > output_regular.dat
		   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output_irregular.dat
		   DEPENDS ${PROJECT_NAME}.main
		   COMMAND $<TARGET_FILE:${PROJECT_NAME}.main> irregular > output_irregular.dat
		   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(${PROJECT_NAME}.plot
		  COMMAND gnuplot ${CMAKE_CURRENT_SOURCE_DIR}/plot.gnu
		  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output_regular.dat ${CMAKE_CURRENT_BINARY_DIR}/output_irregular.dat ${CMAKE_CURRENT_SOURCE_DIR}/plot.gnu
		  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_target_properties(${PROJECT_NAME}.plot PROPERTIES FOLDER "projects/hldo_sphere")

add_dependencies(projects.hldo_sphere.plots ${PROJECT_NAME}.plot)
