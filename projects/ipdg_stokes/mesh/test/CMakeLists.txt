include(GoogleTest)

set(sources mesh_tests.cc)

add_executable(projects.ipdg_stokes.mesh.test ${sources})
target_link_libraries(projects.ipdg_stokes.mesh.test PUBLIC
		      Eigen3::Eigen GTest::gtest_main lf.assemble lf.mesh.hybrid2d lf.mesh lf.mesh.utils lf.geometry lf.base Boost::program_options projects.ipdg_stokes.assemble)
set_target_properties(projects.ipdg_stokes.mesh.test PROPERTIES FOLDER "projects/ipdg_stokes")
target_compile_features(projects.ipdg_stokes.mesh.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.ipdg_stokes.mesh.test)
