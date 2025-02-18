# Quick Reference - Geometry {#quick_reference_geometry}

[TOC]

> [!caution]
> The content of this page is discussed in @lref_link{sec:parmFE}. Please read before using quick reference.

## Overview

LehrFEM++ provides an interface for geometry data: `lf::geometry::Geometry`. Geometry data for entities is stored as a mapping \f$ \Phi_K(\hat{x}) \f$ from the reference element (unit triangle, square, unit segment) to any general element. \f$ \Phi_K^{*}(\hat{x}) \f$ represents the inverse mapping from the general element to the reference element.

For any one element, it is enough to store the mapping \f$ \Phi_K(\hat{x}) \f$ from the reference element and the type of the [reference element](@ref lf::base::RefEl).

![(Affine) Geometry Mapping for Triangles](manim/parametric_fe_geometry.gif)

More details and mathematical background can be found in @lref_link{sec:parmFE}.

## Geometry Interface {#geometry_interface}

To get the geometry of an entity:

```cpp
for (const lf::mesh::Entity* entity : mesh.Entities(codim)) {
    // Get the geometry of an entity
    const lf::geometry::Geometry* geometry = entity->Geometry();
}
```

### Utility Functions {#geometry_utility}

A number of convenience functions are provided by the `Geometry` class.

```cpp
const lf::geometry::Geometry* geometry = entity->Geometry();

// Get the volume of an entity (for edges length, for triangles area, etc.)
double v = lf::geometry::Volume(*geometry);

// Get the corners of an entity as a dxn matrix
// where d is the dimension of physical space, and n the number of corners.
// (for edges 2 corners, for triangles 3 corners, etc.)
Eigen::MatrixXd corners = lf::geometry::Corners(*geometry);
```

In addition, we can test for some properties of the geometry:

```cpp
// Test if geometry is a non-degenerate bilinear quadrilateral (corners matrix with 4 cols)
bool is_not_deg = lf::geometry::assertNonDegenerateQuad(corners);

// Test if geometry is a non-degenerate triangle (corners matrix with 3 cols)
bool is_not_deg = lf::geometry::assertNonDegenerateTriangle(corners);
```

## Geometry Methods {#geometry_methods}

Objects implementing the `lf::geometry::Geometry` interface provide the following methods:

```cpp
// Dimension of the local coordinate system
unsigned dim_local = geometry->DimLocal(); 

// Dimension of the physical coordinate system
unsigned dim_global = geometry->DimGlobal();

// Get the reference element for the geometry
lf::base::RefEl ref_el = geometry->RefEl();

// Check if mapping is affine
bool is_affine = geometry->isAffine();
```

## Transformations and Mappings {#transformations}

This part is discussed in detail in @lref_link{sec:parmFE}. Usage of the following methods is discussed in the [Assembly Quick Reference](@ref quick_reference_assembly).

> [!note] 
> The following methods accept a matrix of points in the local coordinate system as input. The points are stored in a matrix where each column represents a point. The number of rows corresponds to the dimension of the local coordinate system.

### Global {#global}

The `lf::geometry::Geometry` class provides the `lf::geometry::Geometry::Global` method to map points from local to global coordinates. In LehrFEM++, `local` points generally refer to points on the reference element.

```cpp
// Define two point on the reference element (local).
Eigen::MatrixXd local_points(2, 2);
local_points << 0.1, 0.5,
                0.6, 0.2;

// Map the points from the local into the global space.
Eigen::MatrixXd global = geometry->Global(local_points);
```

![Mapping of points from local to global coordinates](manim/mapping_global.gif)

### IntegrationElement {#integration_element}

The method `lf::geometry::Geometry::IntegrationElement` computes the integration element 

\f[
    g(\xi) := \sqrt{\mathrm{det}\left|D\Phi^T(\xi) D\Phi(\xi) \right|}
\f]

for each point in `points` as an `Eigen::VectorXd`.

```cpp
// Compute the integration element for each point in points
Eigen::VectorXd integration_element = geometry->IntegrationElement(local_points);
```

> [!note]
> The following two methods return a matrix for each evaluated point. For efficiency reasons, the returned matrices are horizontally stacked. The number of rows corresponds to the dimension of the global coordinate system.

### Jacobian {#jacobian}

The method `lf::geometry::Geometry::Jacobian` computes the Jacobian matrix \f$ D\Phi(\hat{x}) \f$ for each point in `points`.

```cpp
// Compute the Jacobian matrix for each point in points
Eigen::MatrixXd jacobian = geometry->Jacobian(local_points);
```

The Jacobian evaluated at every point is itself a matrix of size `dim_global x num_points`. The JacobianInverseGramian evaluated at every point is a matrix of size `dim_local x (dim_global * num_points)`. To access the Jacobian at a specific point, use the following code:

```cpp
unsigned dim_global = geometry->DimGlobal();
unsigned dim_local = geometry->DimLocal();
// Access the Jacobian matrix for the n-th point (starting at 0)
Eigen::MatrixXd jacobian_n = jacobian.block(0, n * dim_local, dim_global, dim_local);
```

If `dim_local == dim_global == 2` and we pass three points for evaluation, [Geometry::Jacobian](@ref lf::geometry::Geometry::Jacobian) returns a \f$ 2 \times 6 \f$ matrix:

![Block matrix returned by Geometry::Jacobian](manim/jacobian_block_access.gif)

To access the Jacobian for the second point, we can use [Eigen block access](https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html):

```cpp
// Access 2x2 Jacobian matrix starting at row 0 and column 2
Eigen::MatrixXd jacobian_2 = jacobian.block(0, 2, 2, 2);
```

### JacobianInverseGramian {#jacobian_inverse_gramian}

The method `lf::geometry::Geometry::JacobianInverseGramian` computes the inverse of the JacobianInverseGramian matrix for each point in `points`. Details can be found in @lref_link{rem:jti} and the [method docs](@ref lf::geometry::Geometry::JacobianInverseGramian). Similarly to [Geometry::Jacobian](@ref lf::geometry::Geometry::Jacobian), it also returns a horizontally stacked matrix. Individual matrices can be accessed using `Eigen` block. 

```cpp
// Compute the inverse of the Jacobian matrix for each point in points
Eigen::MatrixXd jacobian_inv = geometry->JacobianInverseGramian(local_points);

// Access the JacobianInverseGramian matrix for the n-th point (starting at 0)
Eigen::MatrixXd jacobian_inv_n = jacobian_inv.block(0, n * dim_local, dim_global, dim_local);
```

<!-- Next and previous buttons -->
<div class="section_buttons">
 
| Previous                          |                                     Next |
| :-------------------------------- | ---------------------------------------: |
| [Mesh](@ref quick_reference_mesh) | [DOFHandlers](@ref quick_reference_dofh) |
 
</div>