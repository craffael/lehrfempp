## LehrFEM++ Parameterization Breakdown Demo

This example code illustrates when the second-order parametrization breaks 
down. More concretely, we investigate when the volume computation becomes 
inexact due to large negative values of the determinant of the jacobian 
evaluated at random points within the geometry element.

The `parametrization_breakdown.cc` file evaluates the `Global()` method and 
computes the determinant of the jacobian at random points for `TriaO2` and 
`QuadO2` objects. To check whether the volume computation is exact, we refine 
all three geometry objects by means of regular refinement and compare the 
refined volume to the original volume. All computations are stored in 
`build/examples/geometry/parametrization_breakdown/results`.

In a second step we visualize the parametrization breakdown in the Jupyter 
Notebook `pararmetrization_breakdown.ipynb`. We plot the random global 
evaluations for the original and refined geometry objects and we also 
visualize the values of the jacobian's determinant on the reference element. 
We observe that for non-degenerate geometry objects the `Global()` evaluations 
lie within the global objects' boundaries and that the jacobian's determinant 
remains positive. However, the degenerate objects fails to conserve their 
volumes throughout refinement and we see that this is caused by large negative 
determinant values. The file `generate_plots.py` stores the plots in 
`build/examples/geometry/parametrization_breakdown/plots` via the command 
`python generate_plots.py /path/to/results`.
