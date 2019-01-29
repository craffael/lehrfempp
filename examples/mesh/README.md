## Test Mesh Output Demo

This demo calls the utility function `lf::mesh::test_utils::GenerateHybrid2DTestMesh()`
implemented in
[`test_meshes.cc`](https://github.com/craffael/lehrfempp/blob/master/lib/lf/mesh/test_utils/test_meshes.cc)
to generate a planar 2D hybrid mesh. This mesh is then output in various formats:

- TikZ output via `lf::mesh::utils::writeTikZ()`. TikZ output has to be processed with
  LaTeX, see `test_meshes.tex`
- Matlab output via `lf::mesh::utils::writeMatlab()`. Use the MATLAB script
  `plot_lf_mesh.m` for generating graphics.
- Pyhton output via `lf::io::writeMatplotlib()`. Visualization can be done using the
  Python script `plot_mesh.py`
- TODO: vtk output via


