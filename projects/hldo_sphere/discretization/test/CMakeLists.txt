include(GoogleTest)

add_executable(projects.hldo_sphere.discretization.test whitney_hodge_laplace_test.cc dirac_operator_test.cc)
target_link_libraries(projects.hldo_sphere.discretization.test PUBLIC
		      Eigen3::Eigen GTest::gtest_main lf.base lf.assemble lf.mesh lf.mesh.hybrid2d lf.mesh.utils lf.geometry projects.hldo_sphere.assemble projects.hldo_sphere.discretization)
set_target_properties(projects.hldo_sphere.discretization.test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(projects.hldo_sphere.discretization.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.hldo_sphere.discretization.test)
