# Donut2d Example

With this example we can demonstrate two things:
1) Approximation of a donut geometry
2) Approximation of a arctan like geometry which shows the superiority 
   of the transfinite approach vs. 2.nd order geometry approximation

In both cases, the code loads the mesh and
- Computes the volume of all mesh elements and compares it to the exact 
  volume of the donut (also for the arctan example)
- Checks that the edges of every mesh element intersect only at the corners.

## 1. Donut geometry
To run this example, just execute the executable with a specific mesh, e.g. 
```
./examples.brep.donut2d --mesh=donut2d_2.45.msh
```
which will run it with a mesh where Gmsh's "Global MeshSize Factor" has been set to 2.45.
(The meshes consist of 2nd order elements)

If we use a larger meshsize factor,
```
./examples.brep.donut2d --mesh=donut2d_2.53.msh
```
We can see that transfinite interpolation as well as the normal second
order elements ("Plain Mesh") yield elements where the edges of some trias
self-intersect. This leads to the volume being calculated wrongly.
However this is not such a big problem because it is easily detected
(edges of the TransfiniteTriangle self intersect), also the edges of the 
original first order elements intersect with the Boundary representation.

--> Here the transfinite interpolation doesn't have an advantage over 2nd order elements.

## 2. ArcTan like Geometry
Command line:
```
./examples.brep.donut2d --mesh=arctan_like.msh --brep=arctan_like.brep
```

If we use an ArcTan-like Geometry with coarse mesh elements, it can happen
that the second-order elements degenerate even though first order elements would not.
This happens mostly because the 2nd-order lagrange basis takes also negative values.

We can see that the transfinite elements master this challenge much better and actually
yields the correct volume of 1.6997 (according to FreeCad).
This is much better because in this case the first order mesh is fine
(no intersection with BRep).

### 2.2 Bonus
It seems that 2nd and 3rd order mesh both don't work for the ArcTan geometry.
In both cases Gmsh outputs a warning about negativ jacobians :)