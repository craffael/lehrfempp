include(GoogleTest)

set(sources
  utils.h
  occt_brep_model_tests.cc
  occt_curve_geometry_tests.cc
)

add_executable(lf.brep.occt.test ${sources})
target_link_libraries(lf.brep.occt.test PUBLIC GTest::gtest_main
  lf.brep.occt
)
set_target_properties(lf.brep.occt.test PROPERTIES FOLDER "lib")
target_compile_features(lf.brep.occt.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.brep.occt.test)