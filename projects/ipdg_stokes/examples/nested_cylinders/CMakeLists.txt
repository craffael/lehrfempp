project(projects.ipdg_stokes.examples.nested_cylinders)

add_executable(${PROJECT_NAME}.main nested_cylinders.cc)
target_link_libraries(${PROJECT_NAME}.main PUBLIC Eigen3::Eigen lf.base lf.mesh lf.mesh.hybrid2d lf.quad lf.refinement lf.io projects.ipdg_stokes.assemble projects.ipdg_stokes.mesh projects.ipdg_stokes.post_processing ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_compile_features(${PROJECT_NAME}.main PUBLIC cxx_std_17)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output_builder.dat
		   DEPENDS $<TARGET_FILE:${PROJECT_NAME}.main>
		   COMMAND ./${PROJECT_NAME}.main builder > output_builder.dat
		   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/output_files.dat
		   DEPENDS $<TARGET_FILE:${PROJECT_NAME}.main>
		   COMMAND ./${PROJECT_NAME}.main files > output_files.dat
		   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(${PROJECT_NAME}.plot
		  COMMAND gnuplot ${CMAKE_CURRENT_SOURCE_DIR}/plot.gnu
		  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/output_builder.dat ${CMAKE_CURRENT_BINARY_DIR}/output_files.dat ${CMAKE_CURRENT_SOURCE_DIR}/plot.gnu
		  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(${PROJECT_NAME}.irregular
		  DEPENDS ${PROJECT_NAME}.main
		  COMMAND ./${PROJECT_NAME}.main irregular > /dev/null
		  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_dependencies(projects.ipdg_stokes.plots ${PROJECT_NAME}.plot ${PROJECT_NAME}.irregular)
