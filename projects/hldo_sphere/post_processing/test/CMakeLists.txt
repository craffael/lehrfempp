include(GoogleTest)

set(sources norms_tests.cc mesh_function_test.cc)

add_executable(projects.hldo_sphere.post_processing.test ${sources})
target_link_libraries(projects.hldo_sphere.post_processing.test PUBLIC
		      Eigen3::Eigen GTest::gtest_main lf.base lf.assemble lf.mesh lf.mesh.hybrid2d lf.mesh.utils lf.uscalfe lf.refinement projects.hldo_sphere.post_processing)
set_target_properties(projects.hldo_sphere.post_processing.test PROPERTIES FOLDER "projects/ipdg_stokes")
target_compile_features(projects.hldo_sphere.post_processing.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.hldo_sphere.post_processing.test)
