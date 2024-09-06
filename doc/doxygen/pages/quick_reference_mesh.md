# Quick Reference - Mesh {#quick_reference_mesh}

[TOC]

Meshes in LehrFEM++ are generally managed through shared pointers. A simple test mesh can be generated using the following code snippet:

```cpp
// Generate a simple test mesh
const std::shared_ptr<const lf ::mesh::Mesh> mesh_ptr = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

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
auto num_entities = mesh_p->NumEntities(0);

// Access number of edges
auto num_nodes = mesh_p->NumEntities(1);

// Get iterable over all cells
auto cells = mesh_p->Entities(0);

// Iterate over all cells
for (const auto* cell : cells) {
    // Do something with the cell
    std::cout << cell->Index() << std::endl;
}
```

<div class="section_buttons">
 
| Previous          |                              Next |
|:------------------|----------------------------------:|
|  | [Geometry](quick_reference_geometry.md) |
 
</div>
