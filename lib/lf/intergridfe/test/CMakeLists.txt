include(GoogleTest)

set(src
    transfer_tests.cc
)

add_executable(lf.intergridfe.tests ${src})
target_link_libraries(lf.intergridfe.tests PUBLIC Eigen3::Eigen GTest::gtest_main lf.mesh lf.mesh.test_utils lf.refinement lf.uscalfe lf.io)
target_compile_features(lf.intergridfe.tests PUBLIC cxx_std_17)
gtest_discover_tests(lf.intergridfe.tests)
