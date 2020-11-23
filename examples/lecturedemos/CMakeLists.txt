# cmake configuration file for compilation of LehrFEM++ lecture demo codes:

add_subdirectory(convergencestudies)

set(lecturedemo
  lecturedemotwonorm.h lecturedemotwonorm.cc
  lecturedemomesh.h lecturedemomesh.cc
  lecturedemodof.h lecturedemodof.cc
  lecturedemoassemble.h lecturedemoassemble.cc
  lecturedemoquad.h lecturedemoquad.cc
  lecturedemorefine.h lecturedemorefine.cc
  lecturedemomeshfunction.h lecturedemomeshfunction.cc
  lecturedemo.h lecturedemomain.cc)

add_executable(examples.lecturedemo ${lecturedemo})

target_link_libraries(examples.lecturedemo
  PUBLIC Eigen3::Eigen
  Boost::program_options
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.io
  lf.quad
  lf.assemble
  )
set_target_properties(examples.lecturedemo PROPERTIES FOLDER "examples")

add_custom_target(examples.lecturedemo_run COMMAND examples.lecturedemo)
add_dependencies(examples_run examples.lecturedemo_run)
set_target_properties(examples.lecturedemo_run PROPERTIES FOLDER "examples")
  

  

