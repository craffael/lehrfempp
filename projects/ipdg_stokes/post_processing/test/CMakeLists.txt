include(GoogleTest)

set(sources norms_tests.cc)

add_executable(projects.ipdg_stokes.post_processing.test ${sources})
target_link_libraries(projects.ipdg_stokes.post_processing.test PUBLIC
		      Eigen3::Eigen GTest::main lf.base lf.assemble lf.mesh lf.mesh.hybrid2d lf.mesh.utils projects.ipdg_stokes.post_processing)
target_compile_features(projects.ipdg_stokes.post_processing.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.ipdg_stokes.post_processing.test)
