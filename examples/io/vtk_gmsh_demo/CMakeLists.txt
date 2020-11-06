# cmake configuration file for compilation of LehrFEM++ examples: io demos

set(vtk_gmsh_demo io_demo.cc)

add_executable(examples.io.vtk_gmsh_demo ${vtk_gmsh_demo})

target_link_libraries(examples.io.vtk_gmsh_demo
  PUBLIC Eigen3::Eigen
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.io
  )
set_target_properties(examples.io.vtk_gmsh_demo PROPERTIES FOLDER "examples")
  
add_custom_target(examples.io.vtk_gmsh_demo_run COMMAND examples.io.vtk_gmsh_demo)
add_dependencies(examples_run examples.io.vtk_gmsh_demo_run)
set_target_properties(examples.io.vtk_gmsh_demo_run PROPERTIES FOLDER "examples")

set(mesh_analysis_demo mesh_analysis_demo.cc)
add_executable(examples.io.mesh_analysis_demo ${mesh_analysis_demo})
target_link_libraries(examples.io.mesh_analysis_demo
  PUBLIC Eigen3::Eigen
  lf.mesh
  lf.mesh.test_utils
  lf.mesh.utils
  lf.mesh.hybrid2d
  lf.refinement
  lf.io
  )
set_target_properties(examples.io.mesh_analysis_demo PROPERTIES FOLDER "examples")

