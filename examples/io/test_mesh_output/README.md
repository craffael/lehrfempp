## LehrFEM++ Mesh Output Demo

The file `test_mesh_output_demo.cc` outputs the topology and geometry of 
pre-defined LehrFEM++ test meshes in various formats:

- **Python format:** as a `*.csv`-file which can be visualized with the 
following command: `python lib/lf/io/plot_mesh.py 
build/examples/io/test_mesh_output/test_mesh_1.csv`.

- **MATLAB format:** as a `*.m`-file. The MATLAB script 
`matlab_testmesh_plot_driver.m` shows how to visualize the mesh. It relies on 
the MATLAB function [`plot_lf_mesh`](
https://github.com/craffael/lehrfempp/blob/master/lib/lf/io/plot_lf_mesh.m
"plot_lf_mesh.m").

- **VTK format:** as a`*.vtk`-file to be viewed with 
[ParaView](https://www.paraview.org/ "paraview.org"). After starting the 
ParaView Client and opening the file, one can view the mesh by clicking `Apply` 
in the `Properties` section and as setting the `Representation` to `Wireframe`.
   
- **(LaTeX) TikZ format:** as a `*.tex`-file which can be included into a 
LaTeX document as illustrated [here](
https://github.com/craffael/lehrfempp/blob/master/examples/io/test_mesh_output/test_meshes.tex
 "LaTeX wrapper file"). You can also view the file directly with the 
[TikzEdt](http://www.tikzedt.org/ "tikzedt.org") editor.
