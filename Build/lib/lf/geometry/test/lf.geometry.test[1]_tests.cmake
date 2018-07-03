add_test( PointTest.checkCoord /u/magina/Documents/lehrfempp/Build/lib/lf/geometry/test/lf.geometry.test [==[--gtest_filter=PointTest.checkCoord]==] --gtest_also_run_disabled_tests)
set_tests_properties( PointTest.checkCoord PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/geometry/test)
set( lf.geometry.test_TESTS PointTest.checkCoord)
