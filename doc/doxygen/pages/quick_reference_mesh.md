# Quick Reference - Mesh {#quick_reference_mesh}

[TOC]

@gh_edit

> [!caution]
> The contents of this page is discussed in @lref_link{sec:meshdata}. Please read this section before using the quick reference.

<!-- Mesh Information and Mesh Data Structures are discussed in detail in @lref_link{sec:meshdata}. -->

## Overview

Meshes in LehrFEM++ are generally managed through shared pointers. A simple test mesh can be generated using the following code snippet:

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Overview

## Mesh Access {#mesh_access}

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Access

The mesh object can be accessed through the pointer `mesh_p` (will be used throughout this page). The mesh object provides access to the following information:

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Access 0

## Entities {#entities}

Entities implement the lf::mesh::Entity interface and are the building blocks of a mesh. LehrFEM++ provides a number of different types of entities:

- **Triangles** lf::mesh::hybrid2d::Triangle
- **Quadrilaterals** lf::mesh::hybrid2d::Quadrilateral
- **Segments** lf::mesh::hybrid2d::Segment
- **Points** lf::mesh::hybrid2d::Point

A complete list can be found in the inheritance diagram of the lf::mesh::Entity interface class. Entities can be accessed by index or by iterating over all entities of a given co-dimension.

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Entities

Entities have a number of properties that can be accessed:

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Entities 0

It is also possible to access sub-entities of an entity, see lf::mesh::Entity::SubEntities for details.

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Entities 1

A slightly more nuanced concept is the relative orientation of sub-entities of the next higher co-dimension. A detailed explanation can be found in the @lref_link{rem:ori}.

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Entities 2

## Mesh Data Sets {#mesh_data_sets}

Mesh data sets are used to store data with entities of the mesh. A fe common use cases are:

- Flags that mark boundary entities.
- Material parameters for mesh elements (codim=0)

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Data Sets

In this example, a mesh data set is created that stores boolean values for each entity/edge of the mesh.

To flag all entities on the boundary of the mesh, the following code snippet can be used:


@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Data Sets 0

See also [Quick Reference - Boundary Conditions](@ref quick_reference_bc) for more details on boundary conditions.

The interface is defined in the lf::mesh::MeshDataSet class. More details can be found in the documentation of the lf::mesh::MeshDataSet class.

## Mesh Functions {#mesh_functions}

Mesh functions are wrappers around a functor that can be evaluated on the entire mesh. The interface is defined in lf::mesh::utils::MeshFunction. A more detailed definition can be found in the @lref_link{supp:mshfn}.

The most general representative of a mesh function is lf::mesh::utils::MeshFunctionGlobal. It takes a functor that can be evaluated on the entire mesh. The functor must provide an operator() which takes a point within the entity reference Element and returns a value. Lambda functions are a common way to define such functors.

For efficiency reasons the evaluation points are passed to mesh functions as the columns of a  \f$d \times n \f$ - matrix (where \f$d\f$ agrees with the local dimension of the entity and \f$n\f$ is the number of points). This allows for the evaluation of multiple points at once. The following code snippet demonstrates how to define a mesh function that evaluates the function

\f[
    \alpha(x) = 0.5 \cdot \|x\|
\f]

on a cell.

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Functions

There are is also a short-hand for mesh functions that are constant everywhere on the mesh: lf::mesh::utils::MeshFunctionConstant.

MeshFunction objects support binary arithmetic operations +,-, and *, including scalar multiplication, provided that such operations are possible for their underlying types. As well as unary operations such as -, transpose(), and squaredNorm().

### Finite Element Mesh Functions {#femf}

A special case of mesh functions are finite element mesh functions. They are used to represent the solution of a finite element problem on a mesh.

<!-- TODO (barmstron): Add note on Integration and evaluation of FE Mesh Functions-->

- lf::fe::MeshFunctionFE
- lf::fe::MeshFunctionGradFE

## Mesh Creation {#mesh_creation}

The standard ways to create a Mesh object are:

1. **Direct generation** via the lf::mesh::MeshFactory::Build() method. (This method is mainly used internally by LehrFEM++). Details can be found in the documentation of the lf::mesh::MeshFactory class. A brief example is outlined below.

2. **Calling LehrFEM++'s generator of test meshes** lf::mesh::test_utils::GenerateHybrid2DTestMesh(). The method provides access to a number of simple predefined test meshes. The following code snippet demonstrates how to generate a simple test mesh:

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Creation

3. **Reading a mesh from file**: see lf::io::GmshReader, and invoking lf::io::GmshReader::mesh(). The following code snippet demonstrates how to read a mesh from a file:

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Creation 0

4. **Refining an existing mesh**, see lf::refinement::MeshHierarchy or the short example below.

### Mesh Refinement {#mesh_refinement}

LehrFEM++ provides a number of mesh refinement tools included in the lf::refinement namespace. 

Mesh refinement using LehrFEM++ is covered in @lref_link{ss:ref} and heavily used in @lref_link{cha:cvg}.

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Refinement

### Mesh Builder {#mesh_builder}

Meshes can be built 'manually' using the lf::mesh::MeshFactory class. It follows a [builder pattern](https://refactoring.guru/design-patterns/builder) and allows the user to add entities to the mesh. 

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Builder

We can also use the mesh builder to create a tensor product mesh. LehrFEM++ provides a number of mesh builders for different types of meshes.

- **lf::mesh::utils::TPTriagMeshBuilder** for a tensor product mesh of triangles.
- **lf::mesh::utils::TPQuadMeshBuilder** for a tensor product mesh of quadrilaterals.
- **lf::mesh::utils::TorusMeshBuilder** for a mesh of a torus.

The following code snippet demonstrates how to create a 100x100 tensor product mesh of triangles:

@snippet{trimleft} quick_reference/qr_mesh_snippets.cpp mesh Mesh Builder 0


<!-- Next and previous buttons -->
<div class="section_buttons">

| Previous |                                      Next |
| :------- | ----------------------------------------: |
|          | [Geometry](@ref quick_reference_geometry) |

</div>
