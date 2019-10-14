include(GoogleTest)

add_executable(projects.ipdg_stokes.assemble.test matrix_assembly_tests.cc vector_assembly_tests.cc)
target_link_libraries(projects.ipdg_stokes.assemble.test PUBLIC
		      Eigen3::Eigen GTest::main lf.base lf.assemble lf.mesh lf.mesh.hybrid2d lf.mesh.utils lf.geometry projects.ipdg_stokes.assemble)
target_compile_features(projects.ipdg_stokes.assemble.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.ipdg_stokes.assemble.test)
