# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(reo_sources ref_el_output.cc)

add_executable(lf.examples.base.ref_el_output ${reo_sources})

target_link_libraries(lf.examples.base.ref_el_output PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.hybrid2d lf.mesh.test_utils lf.mesh.utils)

target_compile_features(lf.examples.base.ref_el_output PUBLIC cxx_std_17)


