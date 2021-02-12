include(GoogleTest)

set(sources
  brep_lagr_segment_tests.cc
  brep_mesh_factory_transfinite_tests.cc
  brep_segment_tests.cc
  brep_surface_segment_tests.cc
  brep_tria_transfinite_tests.cc
  curve_straight_segment_tests.cc
  utils.h
)

add_executable(lf.brep.geom.test ${sources})
target_link_libraries(lf.brep.geom.test PUBLIC GTest::gtest_main GTest::gmock_main
  lf.brep.geom
  lf.brep.test_utils
  lf.geometry
  lf.geometry.test_utils
  lf.quad
)
set_target_properties(lf.brep.geom.test PROPERTIES FOLDER "lib")
target_compile_features(lf.brep.geom.test PUBLIC cxx_std_17)
gtest_discover_tests(lf.brep.geom.test)