include(GoogleTest)

set(sources
  curve_circle_tests.cc
)

add_executable(lf.brep.test_utils.test ${sources})
target_link_libraries(lf.brep.test_utils.test PUBLIC GTest::gtest_main
  lf.brep.test_utils
)
set_target_properties(lf.brep.test_utils.test PROPERTIES FOLDER "lib")
target_compile_features(lf.brep.test_utils.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.brep.test_utils.test)