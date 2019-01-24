# cmake configuration file for compilation of LehrFEM++ lecture demo codes:

set(lecturedemo lecturedemomesh.h lecturedemomesh.cc lecturedemomain.cc)

add_executable(examples.lecturedemo ${lecturedemo})

target_link_libraries(examples.lecturedemo
  PUBLIC Eigen3::Eigen Boost::filesystem Boost::system
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.io
  )
  
    

  

