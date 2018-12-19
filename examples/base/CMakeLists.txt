# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(reo_sources ref_el_output.cc)

add_executable(examples.base.ref_el_output ${reo_sources})

target_link_libraries(examples.base.ref_el_output PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.test_utils lf.mesh.utils)

set(sv_demo_sources static_vars_demo.cc)

add_executable(examples.base.static_vars_demo ${sv_demo_sources})

target_link_libraries(examples.base.static_vars_demo PUBLIC Eigen3::Eigen Boost::boost GTest::main lf.mesh.test_utils lf.mesh.utils lf.geometry lf.refinement)

set(bpo_demo_sources boost_po_demo.cc)
add_executable(lf.examples.base.boost_po_demo ${bpo_demo_sources})
target_link_libraries(lf.examples.base.boost_po_demo PUBLIC
                      Boost::program_options
                      lf.mesh.utils)

add_subdirectory(comm)
