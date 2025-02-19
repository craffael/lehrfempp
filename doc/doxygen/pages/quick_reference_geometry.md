# Quick Reference - Geometry {#quick_reference_geometry}

[TOC]

@gh_edit

> [!caution]
> The content of this page is discussed in @lref_link{sec:parmFE}. Please read before using quick reference.

## Overview

LehrFEM++ provides an interface for geometry data: `lf::geometry::Geometry`. Geometry data for entities is stored as a mapping \f$ \Phi_K(\hat{x}) \f$ from the reference element (unit triangle, square, unit segment) to any general element. \f$ \Phi_K^{*}(\hat{x}) \f$ represents the inverse mapping from the general element to the reference element.

For any one element, it is enough to store the mapping \f$ \Phi_K(\hat{x}) \f$ from the reference element and the type of the [reference element](@ref lf::base::RefEl).

![(Affine) Geometry Mapping for Triangles](manim/parametric_fe_geometry.gif)

More details and mathematical background can be found in @lref_link{sec:parmFE}.

## Geometry Interface {#geometry_interface}

To get the geometry of an entity:

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Geometry Interface

### Utility Functions {#geometry_utility}

A number of convenience functions are provided by the `Geometry` class.

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Utility Functions

In addition, we can test for some properties of the geometry:

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Utility Functions 0

## Geometry Methods {#geometry_methods}

Objects implementing the `lf::geometry::Geometry` interface provide the following methods:

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Geometry Methods

## Transformations and Mappings {#transformations}

This part is discussed in detail in @lref_link{sec:parmFE}. Usage of the following methods is discussed in the [Assembly Quick Reference](@ref quick_reference_assembly).

> [!note] 
> The following methods accept a matrix of points in the local coordinate system as input. The points are stored in a matrix where each column represents a point. The number of rows corresponds to the dimension of the local coordinate system.

### Global {#global}

The `lf::geometry::Geometry` class provides the `lf::geometry::Geometry::Global` method to map points from local to global coordinates. In LehrFEM++, `local` points generally refer to points on the reference element.

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Global

![Mapping of points from local to global coordinates](manim/mapping_global.gif)

### IntegrationElement {#integration_element}

The method `lf::geometry::Geometry::IntegrationElement` computes the integration element 

\f[
    g(\xi) := \sqrt{\mathrm{det}\left|D\Phi^T(\xi) D\Phi(\xi) \right|}
\f]

for each point in `points` as an `Eigen::VectorXd`.

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry IntegrationElement

> [!note]
> The following two methods return a matrix for each evaluated point. For efficiency reasons, the returned matrices are horizontally stacked. The number of rows corresponds to the dimension of the global coordinate system.

### Jacobian {#jacobian}

The method `lf::geometry::Geometry::Jacobian` computes the Jacobian matrix \f$ D\Phi(\hat{x}) \f$ for each point in `points`.

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Jacobian

The Jacobian evaluated at every point is itself a matrix of size `dim_global x num_points`. The JacobianInverseGramian evaluated at every point is a matrix of size `dim_local x (dim_global * num_points)`. To access the Jacobian at a specific point, use the following code:

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Jacobian 0

If `dim_local == dim_global == 2` and we pass three points for evaluation, [Geometry::Jacobian](@ref lf::geometry::Geometry::Jacobian) returns a \f$ 2 \times 6 \f$ matrix:

![Block matrix returned by Geometry::Jacobian](manim/jacobian_block_access.gif)

To access the Jacobian for the second point, we can use [Eigen block access](https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html):

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry Jacobian 1

### JacobianInverseGramian {#jacobian_inverse_gramian}

The method `lf::geometry::Geometry::JacobianInverseGramian` computes the inverse of the JacobianInverseGramian matrix for each point in `points`. Details can be found in @lref_link{rem:jti} and the [method docs](@ref lf::geometry::Geometry::JacobianInverseGramian). Similarly to [Geometry::Jacobian](@ref lf::geometry::Geometry::Jacobian), it also returns a horizontally stacked matrix. Individual matrices can be accessed using `Eigen` block. 

@snippet{trimleft} quick_reference/qr_geometry_snippets.cpp geometry JacobianInverseGramian

<!-- Next and previous buttons -->
<div class="section_buttons">
 
| Previous                          |                                     Next |
| :-------------------------------- | ---------------------------------------: |
| [Mesh](@ref quick_reference_mesh) | [DOFHandlers](@ref quick_reference_dofh) |
 
</div>