add_test( lf_mesh.buildStructuredMesh /u/magina/Documents/lehrfempp/Build/lib/lf/mesh/test/lf.mesh.test [==[--gtest_filter=lf_mesh.buildStructuredMesh]==] --gtest_also_run_disabled_tests)
set_tests_properties( lf_mesh.buildStructuredMesh PROPERTIES WORKING_DIRECTORY /u/magina/Documents/lehrfempp/Build/lib/lf/mesh/test)
set( lf.mesh.test_TESTS lf_mesh.buildStructuredMesh)
