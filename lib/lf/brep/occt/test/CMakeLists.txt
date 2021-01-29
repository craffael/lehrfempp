include(GoogleTest)

set(sources
  occt_brep_curve_tests.cc
  occt_brep_model_tests.cc
  occt_brep_surface_tests.cc
  utils.h
)

add_executable(lf.brep.occt.test ${sources})
target_link_libraries(lf.brep.occt.test PUBLIC GTest::gtest_main
  lf.brep.occt
  lf.brep.test_utils
)
set_target_properties(lf.brep.occt.test PROPERTIES FOLDER "lib")
target_compile_features(lf.brep.occt.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.brep.occt.test)