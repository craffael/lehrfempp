## LehrFEM++ mesh output demo

The small example code in this folder outputs the topology and geometry of pre-defined
LehrFEM++ test meshes in varous formats:

- In "pyhton format" as a `*.py`-file to be parsed by a special Python function:

   TODO: describe

- In MATLAB format, as a `*.m`-file. The MATLAB script `matlab_testmesh_plot_driver.m`
shows how to visualize the mesh. It relies on the MATLAB function `plot_lf_mesh`
implemented
[here](https://github.com/craffael/lehrfempp/blob/master/lib/lf/io/plot_lf_mesh.m
"plot_lf_mesh.m").

- VTK output, stored in an `*.vtk`-file. This has to be postprocessed with Paraview. 

   TODO: Describe use of Paraview
   
- LaTeX-Tikz output as a `*.tex`-file, which can be included into a LaTeX docuemnt, see 
  [the
  file](https://github.com/craffael/lehrfempp/blob/master/examples/io/test_mesh_output/test_meshes.tex
  "LaTeX wrapper file") for a sample LaTeX wrapper file.
