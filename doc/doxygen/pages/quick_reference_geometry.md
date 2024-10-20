# Quick Reference - Geometry {#quick_reference_geometry}

[TOC]

LehrFEM++ provides an Interface for geometry data: lf::geometry::Geometry. Geometry data is stored as a mapping from the reference element (unit triangle, square) to the physical element. Details and mathematical background can be found in the [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf) @lref{sec:parmFE}.

## Geometry Interface {#geometry_interface}

To get the geometry of an entity :

```cpp
for (const lf::mesh::Entity* entity : mesh.Entities(codim)) {
    // Get the geometry of an entity
    auto geometry = entity->Geometry();
}
```

### Utility Functions {#geometry_utility}

A number of convenience functions are provided by the Geometry class.

```cpp
auto geometry = entity->Geometry();

// Get the Volume of an entity
double v = lf::geometry::Volume(*geometry);

// Get the corners of an entity
auto corners = lf::geometry::Corners(*geometry);

// Geometry mappings can be composed (This feature is not yet implemented)
// auto geometry2 = entity2->Geometry();
// auto composed_geo = lf::geometry::Compose(geometry, geometry2);
```

In addition we can test for some properties of the geometry:

```cpp
// Get corner coordinates of the entity as a dxn matrix
// where d is the dimension of physical space, and n the number of corners.
auto corners = lf::geometry::Corners(*geometry);

// Test if geometry is a non-degenerate bilinear quadrilateral (corners matrix with 4 cols)
bool is_not_deg = lf::geometry::assertNonDegenerateQuad(corners);

// Test if geometry is a non-degenerate triangle (corners matrix with 3 cols)
bool is_not_deg = lf::geometry::assertNonDegenerateTriangle(corners);
```

## Geometry Methods {#geometry_methods}

Objects implementing the lf::geometry::Geometry interface provide the following methods:

```cpp
// Dimension of the local coordinate system
auto dim_local = geometry->DimLocal(); 

// Dimension of the physical coordinate system
auto dim_global = geometry->DimGlobal();

// Get the reference element for the geometry
auto ref_el = geometry->RefEl();

// Check if mapping is affine
bool is_affine = geometry->IsAffine();
```

## Transformations and Mappings {#transformations}

This part is discussed in detail in [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf) @lref{sec:parmFE}. The lf::geometry::Geometry class provides the following methods to compute and use the underlying mappings:

```cpp
// A number of point in the reference element
auto points = Eigen::MatrixXd(2, 3);
points << 0.0, 1.0, 0.0,
          0.0, 0.0, 1.0;

// Map a the points from the local into the global coordinate system.
auto eval = geometry->Global(points);

// TODO: Jacobian, JacobianInverseGramian, IntegrationElement  

```

<!-- Next and previous buttons -->
<div class="section_buttons">
 
| Previous                        |                            Next |
| :------------------------------ | ------------------------------: |
| [Mesh](quick_reference_mesh.md) | [Dofs](quick_reference_dofs.md) |
 
</div>
