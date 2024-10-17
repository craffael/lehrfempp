# Quick Reference - Mesh {#quick_reference_mesh}

[TOC]

Meshes in LehrFEM++ are generally managed through shared pointers. A simple test mesh can be generated using the following code snippet:

```cpp
// Generate a simple test mesh
const std::shared_ptr<const lf::mesh::Mesh> mesh_ptr =
    lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

// Auto can be used to simplify the type declaration
auto mesh_ptr = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
```

## Mesh Creation

The standard ways to create a Mesh object are:

1. **Direct generation** via the lf::mesh::MeshFactory::Build() method. (This method is mainly used internally by LehrFEM++). Details can be found in the documentation of the lf::mesh::MeshFactory class.

2. **Calling LehrFEM++'s generator of test meshes** lf::mesh::test_utils::GenerateHybrid2DTestMesh(). The method provides access to a number of simple predefined test meshes. The following code snippet demonstrates how to generate a simple test mesh:

    ```cpp
    // Generate a simple test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
    ```

3. **Reading a mesh from file**: see lf::io::GmshReader, and invoking lf::io::GmshReader::mesh(). The following code snippet demonstrates how to read a mesh from a file:

    ```cpp
    // Read a mesh from a file
    auto mesh_p = lf::io::GmshReader().mesh("path/to/mesh.msh");
    ```

4. **Refining an existing mesh**, see lf::refinement::MeshHierarchy.

## Mesh Access

```cpp
// create a pointer to the mesh
auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
```

The mesh object can be accessed through the pointer `mesh_p`. The mesh object provides access to the following information:

```cpp
// Access number of cells
auto num_entities = mesh->NumEntities(0);

// Access number of edges
auto num_nodes = mesh->NumEntities(1);

// Iterate over all cells
for (const auto* cell : mesh->Entities(0)) {
// Do something with entity e.g.: print index of cell
std::cout << mesh->Index(*cell) << std::endl;
}
```

## Mesh Utils

LehrFEM++ provides a number of utility functions to work with meshes.

### Mesh Data Sets 

Mesh data sets are used to store data with entities of the mesh. A common use cases are:

- Flags that mark boundary entities.
- Material parameters for mesh elements (codim=0)

```cpp
// Create a mesh data set storing boolean values for each entity
auto mesh_data_set = lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p);

// Create a MDS storing boolean values for each entity of co-dimension 1 (edges)
auto mesh_data_set_edges = lf::mesh::utils::CodimMeshDataSet<bool>(mesh_p, 1);

// Mesh data sets can be initialized with a default value
auto mesh_data_set = lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p, false);

// Access a (modifiable) value associated with an entity
auto entity = mesh_p->EntityByIndex(0, 0);
auto value = mesh_data_set(*entity);
```

In this example, a mesh data set is created that stores boolean values for each entity/edge of the mesh.


To flag all entities on the boundary of the mesh, the following code snippet can be used:
```cpp
// Create MeshDataSet storing boolean values for each entity indicating if the entity is on the boundary
auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p)};
```cpp


The interface is defined in the lf::mesh::MeshDataSet class. More details can be found in the documentation of the lf::mesh::MeshDataSet class.


<!-- Next and previous buttons -->
<div class="section_buttons">

| Previous          |                              Next |
|:------------------|----------------------------------:|
|  | [Geometry](quick_reference_geometry.md) |

</div>
````
