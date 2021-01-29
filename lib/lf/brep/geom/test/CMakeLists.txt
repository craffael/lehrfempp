include(GoogleTest)

set(sources
  brep_lagr_segment_tests.cc
)

add_executable(lf.brep.geom.test ${sources})
target_link_libraries(lf.brep.geom.test PUBLIC GTest::gtest_main
  lf.brep.geom
  lf.brep.test_utils
  lf.geometry
  lf.quad
)
set_target_properties(lf.brep.geom.test PROPERTIES FOLDER "lib")
target_compile_features(lf.brep.geom.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.brep.geom.test)