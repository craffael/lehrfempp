# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(reo_sources ref_el_output.cc)

add_executable(examples.base.ref_el_output ${reo_sources})

target_link_libraries(examples.base.ref_el_output PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.test_utils lf.mesh.utils)

add_subdirectory(comm)
