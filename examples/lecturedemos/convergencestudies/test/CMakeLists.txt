include(GoogleTest)

add_executable(examples.lecturedemos.convergencestudies.test fespacehp_tests.cc)
target_link_libraries(examples.lecturedemos.convergencestudies.test PUBLIC
    Eigen3::Eigen
    GTest::gtest_main
    lf.base
    lf.assemble
    lf.mesh
    lf.mesh.hybrid2d
    lf.mesh.utils
    lf.uscalfe
    lf.refinement
    lf.io
)
target_compile_features(examples.lecturedemos.convergencestudies.test PUBLIC cxx_std_17)
gtest_discover_tests(examples.lecturedemos.convergencestudies.test)
