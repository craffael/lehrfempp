set(parametrization_breakdown_sources parametrization_breakdown.cc)
add_executable(
        examples.geometry.parametrization_breakdown ${parametrization_breakdown_sources}
)
target_link_libraries(examples.geometry.parametrization_breakdown PUBLIC
        Eigen3::Eigen
        lf.geometry lf.quad lf.refinement lf.mesh.hybrid2d
)
set_target_properties(examples.geometry.parametrization_breakdown PROPERTIES FOLDER "examples")

add_custom_target(examples.geometry.parametrization_breakdown_run COMMAND examples.geometry.parametrization_breakdown)
add_dependencies(examples_run examples.geometry.parametrization_breakdown_run)
set_target_properties(examples.geometry.parametrization_breakdown_run PROPERTIES FOLDER "examples")
