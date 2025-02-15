# Quick Reference - Mesh {#quick_reference_mesh}

[TOC]

> [!caution]
> The contents of this page is discussed in @lref_link{sec:meshdata}. Please read before using quick reference.

<!-- Mesh Information and Mesh Data Structures are discussed in detail in @lref_link{sec:meshdata}. -->

## Overview

Meshes in LehrFEM++ are generally managed through shared pointers. A simple test mesh can be generated using the following code snippet:

```cpp
// Generate a simple test mesh
const std::shared_ptr<const lf::mesh::Mesh> mesh_ptr =
    lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

// Auto can be used to simplify the type declaration
auto mesh_ptr = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
```

## Mesh Access {#mesh_access}

```cpp
// create a pointer to the mesh
const std::shared_ptr<const lf::mesh::Mesh> mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
```

The mesh object can be accessed through the pointer `mesh_p` (will be used throughout this page). The mesh object provides access to the following information:

```cpp
// Access number of cells
unsigned num_entities = mesh_p->NumEntities(0);

// Access number of edges
unsigned num_nodes = mesh_p->NumEntities(1);

```

## Entities {#entities}

Entities implement the lf::mesh::Entity interface and are the building blocks of a mesh. LehrFEM++ provides a number of different types of entities:

- **Triangles** lf::mesh::hybrid2d::Triangle
- **Quadrilaterals** lf::mesh::hybrid2d::Quadrilateral
- **Segments** lf::mesh::hybrid2d::Segment
- **Points** lf::mesh::hybrid2d::Point

A complete list can be found in the inheritance diagram of the lf::mesh::Entity interface class. Entities can be accessed by index or by iterating over all entities of a given co-dimension.

```cpp
// Get entity pointer by index from a mesh
const lf::mesh::Entity* entity = mesh_p->EntityByIndex(0, 0);

// Iterate over all entities of co-dimension 0 (cells)
for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    // Do something with entity e.g.: print index of cell
    std::cout << mesh_p->Index(*cell) << std::endl;
}
```

Entities have a number of properties that can be accessed:

```cpp
// Get the co-dimension of the entity
unsigned codim = entity->Codim();

// Get the geometry of the entity, more details can be found 
// in the Geometry Quick reference
const lf::geometry::Geometry* geometry = entity->Geometry();

// Get the reference element of the entity (e.g. the unit triangle for a general triangle)
lf::base::RefEl ref_el = entity->RefEl();
```

It is also possible to access sub-entities of an entity, see lf::mesh::Entity::SubEntities for details.

```cpp
// Return all sub entities of this entity that have the given 
// co-dimension (w.r.t. this entity!)
// For example, for a cell, the sub-entities of co-dimension 1
std::span<const lf::mesh::Entity* const> sub_entities =
      entity->SubEntities(1);

// If you require a vector of sub-entities, you can use the following code snippet
std::vector<const lf::mesh::Entity*> sub_entities_vec{
    entity->SubEntities(1).begin(), entity->SubEntities(1).end()};

```

A slightly more nuanced concept is the relative orientation of sub-entities of the next higher co-dimension. A detailed explanation can be found in the @lref_link{rem:ori}.

```cpp
// return span of relative orientations of sub-entities of the next higher co-dimension.
std::span<const lf::mesh::Orientation> orientations =
      entity->RelativeOrientations();

// If you require a vector of orientations, you can use the following code snippet
std::vector<lf::mesh::Orientation> orientations_vec{
      entity->RelativeOrientations().begin(),
      entity->RelativeOrientations().end()};
```

## Mesh Data Sets {#mesh_data_sets}

Mesh data sets are used to store data with entities of the mesh. A common use cases are:

- Flags that mark boundary entities.
- Material parameters for mesh elements (codim=0)

```cpp
// Create a mesh data set storing boolean values for each entity
lf::mesh::utils::AllCodimMeshDataSet<bool> mesh_data_set =
      lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p);

// Create a MDS storing boolean values for each entity of co-dimension 1 (edges)
lf::mesh::utils::CodimMeshDataSet<bool> mesh_data_set_edges =
      lf::mesh::utils::CodimMeshDataSet<bool>(mesh_p, 1);

// Mesh data sets can be initialized with a default value
lf::mesh::utils::AllCodimMeshDataSet<bool> mesh_data_set_2 =
      lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p, false);

// Access a (modifiable) value associated with an entity
const lf::mesh::Entity* entity = mesh_p->EntityByIndex(0, 0);
bool value = mesh_data_set(*entity);
```

In this example, a mesh data set is created that stores boolean values for each entity/edge of the mesh.

To flag all entities on the boundary of the mesh, the following code snippet can be used:


```cpp
// Create MeshDataSet storing boolean values for each entity indicating if the entity is on the boundary
lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
```

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

```cpp
// Define a lambda function that takes a point in the reference element of an entity and returns a value
auto alpha = [](Eigen::Vector2d x) -> double {
    return 0.5 * x.norm();
};

// Wrap the lambda function in a MeshFunctionGlobal object
// Type of mesh function template parameter is automatically deduced
lf::mesh::utils::MeshFunctionGlobal mesh_function = lf::mesh::utils::MeshFunctionGlobal(alpha);

// Evaluate the mesh function on a cell
const lf::mesh::Entity* cell = mesh_p->EntityByIndex(0, 0);
std::vector<double> values = mesh_function(*cell, Eigen::Vector2d{0.5, 0.5});

// We only asked for one point, so the vector has only one entry
std::cout << "Value: " << values[0] << std::endl; // Value: 0.790569
```

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

    ```cpp
    // Generate a simple test mesh
    std::shared_ptr<lf::mesh::Mesh> mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
    ```

3. **Reading a mesh from file**: see lf::io::GmshReader, and invoking lf::io::GmshReader::mesh(). The following code snippet demonstrates how to read a mesh from a file:

    ```cpp
    // Read a mesh from a file
    std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader gmsh_reader =
        lf::io::GmshReader(std::move(mesh_factory), "triangle.msh");
    
    // Get pointer to the mesh
    std::shared_ptr<lf::mesh::Mesh> mesh = gmsh_reader.mesh();
    ```

4. **Refining an existing mesh**, see lf::refinement::MeshHierarchy or the short example below.

### Mesh Refinement {#mesh_refinement}

LehrFEM++ provides a number of mesh refinement tools included in the lf::refinement namespace. 

Mesh refinement using LehrFEM++ is covered in @lref_link{ss:ref} and heavily used in @lref_link{cha:cvg}.

```cpp
std::shared_ptr<lf::mesh::Mesh> mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
// Generate mesh hierarchy by uniform refinement with 6 levels

std::shared_ptr<lf::refinement::MeshHierarchy> mesh_seq_p{
    lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, 3)};

// Access the mesh at level 3
std::shared_ptr<lf::mesh::Mesh> mesh_level_3 = mesh_seq_p->getMesh(3);

// We can refine a mesh further by calling
mesh_seq_p->RefineRegular();

// Access the mesh at level 4
std::shared_ptr<lf::mesh::Mesh> mesh_level_4 = mesh_seq_p->getMesh(4);
```

### Mesh Builder {#mesh_builder}

Meshes can be built 'manually' using the lf::mesh::MeshFactory class. It follows a [builder pattern](https://refactoring.guru/design-patterns/builder) and allows the user to add entities to the mesh. 

```cpp
// builder for a hybrid mesh in a world of dimension 2
std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory =
    std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

// Add points (needs to be done before adding entities)
mesh_factory->AddPoint(Eigen::Vector2d{0, 0});    // (0)
mesh_factory->AddPoint(Eigen::Vector2d{1, 0});    // (1)
mesh_factory->AddPoint(Eigen::Vector2d{1, 1});    // (2)
mesh_factory->AddPoint(Eigen::Vector2d{0, 1});    // (3)
mesh_factory->AddPoint(Eigen::Vector2d{0.5, 1});  // (4)

// Add a triangle
// First set the coordinates of its nodes:
Eigen::MatrixXd nodesOfTria(2, 3);
nodesOfTria << 1, 1, 0.5, 0, 1, 1;
mesh_factory->AddEntity(
    lf::base::RefEl::kTria(),  // we want a triangle
    std::array<lf::mesh::Mesh::size_type, 3>{
        {1, 2, 4}},  // indices of the nodes
    std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

// Add a quadrilateral
Eigen::MatrixXd nodesOfQuad(2, 4);
nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
mesh_factory->AddEntity(
    lf::base::RefEl::kQuad(), // we want a quadrilateral
    std::array<lf::mesh::Mesh::size_type, 4>{
        {0, 1, 4, 3}}, // indices of the nodes
    std::make_unique<lf::geometry::QuadO1>(nodesOfQuad)); // node coords

// Build the mesh
std::shared_ptr<lf::mesh::Mesh> mesh = mesh_factory->Build();
```

We can also use the mesh builder to create a tensor product mesh. LehrFEM++ provides a number of mesh builders for different types of meshes.

- **lf::mesh::utils::TPTriagMeshBuilder** for a tensor product mesh of triangles.
- **lf::mesh::utils::TPQuadMeshBuilder** for a tensor product mesh of quadrilaterals.
- **lf::mesh::utils::TorusMeshBuilder** for a mesh of a torus.

The following code snippet demonstrates how to create a 100x100 tensor product mesh of triangles:

```cpp
  // Create a LehrFEM++ square tensor product mesh
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));

  // Set mesh parameters following the builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{99, 99})
      .setNumXCells(100)
      .setNumYCells(100);

  std::shared_ptr<lf::mesh::Mesh> mesh_p = builder.Build();
```


<!-- Next and previous buttons -->
<div class="section_buttons">

| Previous |                                      Next |
| :------- | ----------------------------------------: |
|          | [Geometry](@ref quick_reference_geometry) |

</div>
