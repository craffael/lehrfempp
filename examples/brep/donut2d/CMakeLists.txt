# cmake configuration file for compilation of LehrFEM++ examples: RefEl

if(LF_BUILD_OCCT)
  set(src brep_donut2d.cc)

  add_executable(examples.brep.donut2d ${src})

  target_link_libraries(examples.brep.donut2d PUBLIC Eigen3::Eigen Boost::boost GTest::gtest_main lf.brep.occt lf.brep.geom lf.io lf.mesh.hybrid2d lf.refinement Boost::program_options)
  set_target_properties(examples.brep.donut2d PROPERTIES FOLDER "examples")

  add_custom_target(examples.brep.donut2d_run COMMAND examples.base.ref_el_output)
  add_dependencies(examples_run examples.brep.donut2d_run)
  set_target_properties(examples.brep.donut2d_run PROPERTIES FOLDER "examples")
endif()