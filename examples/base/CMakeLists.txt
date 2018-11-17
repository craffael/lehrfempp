# cmake configuration file for compilation of LehrFEM++ examples: RefEl

set(reo_sources ref_el_output.cc)

add_executable(examples.base.ref_el_output ${reo_sources})

target_link_libraries(examples.base.ref_el_output PUBLIC Eigen3::Eigen Boost::boost GTest::main mesh.test_utils mesh.utils)

set(sv_demo_sources static_vars_demo.cc)

add_executable(examples.base.static_vars_demo ${sv_demo_sources})

target_link_libraries(examples.base.static_vars_demo PUBLIC Eigen3::Eigen Boost::boost GTest::main mesh.test_utils mesh.utils geometry refinement)

  
    

  

