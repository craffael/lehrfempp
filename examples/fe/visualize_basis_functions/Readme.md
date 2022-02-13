This example creates a `out.vtk` file which visualizes all basis functions of a `lf::fe::ScalarFESpace` using Paraview high-order cell types.

**Note** To get around a [bug in Paraview](https://gitlab.kitware.com/paraview/paraview/-/issues/20837)
the first vtk-dataset written to `out.vtk` is a complex trigonometric function so that the tesselate filter works as expected.
This dataset has the name `test`.

All the other vtk-datasets have the names `0`, `1`, ..., `N-1` where `N` is the number of basis functions in the `ScalarFESpace`.

Tips to visualize the file in Paraview:
1) To properly visualize the higher-order elements in Paraview you need to apply the tesselate filter and set a field
   error of e.g. `0.00001` and possibly increase the maximum allowed number of subdivisions.
   See [the documentation of VtkWriter](https://craffael.github.io/lehrfempp/classlf_1_1io_1_1_vtk_writer.html) for more information about the tesselate filter.

2) If you want to extrude the triangle in the z-direction according to the values it takes, you can use the "Warp by Scalar" filter of Paraview.

