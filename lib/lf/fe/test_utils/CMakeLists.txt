
set(sources
  complex_fe_space.h
  complex_fe_space.cc
  test_utils.h
)

lf_add_library(lf.fe.test_utils ${sources})
target_link_libraries(lf.fe.test_utils PUBLIC GTest::gtest_main lf.fe lf.uscalfe)
set_target_properties(lf.fe.test_utils PROPERTIES FOLDER "lib")
