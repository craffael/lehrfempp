set(parametrization_breakdown_sources parametrization_breakdown.cc)
add_executable(
        examples.geometry.parametrization_breakdown ${parametrization_breakdown_sources}
)
target_link_libraries(examples.geometry.parametrization_breakdown PUBLIC
        Eigen3::Eigen
        lf.geometry lf.quad lf.refinement
        stdc++fs
)
target_compile_features(examples.geometry.parametrization_breakdown PUBLIC cxx_std_17)
