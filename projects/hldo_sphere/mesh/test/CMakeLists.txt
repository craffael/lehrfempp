include(GoogleTest)

set(sources mesh_tests.cc)

add_executable(projects.hldo_sphere.mesh.test ${sources})
target_link_libraries(projects.hldo_sphere.mesh.test PUBLIC
		      Eigen3::Eigen GTest::gtest_main lf.assemble lf.mesh.hybrid2d lf.mesh lf.mesh.utils lf.geometry lf.base Boost::program_options projects.hldo_sphere.mesh)
set_target_properties(projects.hldo_sphere.mesh.test PROPERTIES FOLDER "projects/hldo_sphere")
target_compile_features(projects.hldo_sphere.mesh.test PUBLIC cxx_std_17)
gtest_discover_tests(projects.hldo_sphere.mesh.test)
