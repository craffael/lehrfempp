include(GoogleTest)

add_executable(projects.hldo_sphere.assemble.test vector_assemby_tests.cc matrix_assembly_two_form_tests.cc matrix_assembly_one_form_tests.cc)
target_link_libraries(projects.hldo_sphere.assemble.test PUBLIC
		      Eigen3::Eigen GTest::gtest_main lf.base lf.assemble lf.mesh lf.mesh.hybrid2d lf.mesh.utils lf.geometry projects.hldo_sphere.assemble)
set_target_properties(projects.hldo_sphere.assemble.test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(projects.hldo_sphere.assemble.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.hldo_sphere.assemble.test)
