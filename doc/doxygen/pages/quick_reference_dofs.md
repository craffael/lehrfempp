# Quick Reference - Degrees of Freedom (DOF) Handlers and Indexing {#quick_reference_dofs}

[TOC]

> [!caution]
> The contents of this page is discussed in @lref_link{sec:parmFE}. Please read before using quick reference.

## Overview
In finite element methods (FEM), handling degrees of freedom (DOFs) is critical for constructing systems of equations. DOF handlers are responsible for assigning local basis functions to global basis functions, which allows for the assembly of global matrices and vectors from element contributions.

### DOF Handler
A DOF handler manages the relationship between the finite element space and the underlying mesh by mapping local shape functions to global shape functions. In **LehrFEM++**, this is done through objects of the type `lf::assemble::DofHandler`.

The key methods provided by the `DofHandler` interface are:
- `NumDofs()`: Returns the total number of global DOFs.
- `NumLocalDofs(const lf::mesh::Entity &)`: Provides the number of DOFs associated with a particular geometric entity.
- `GlobalDofIndices(const lf::mesh::Entity &)`: Returns the global indices of DOFs associated with a given entity.
- `NumInteriorDofs(const lf::mesh::Entity &)`: Specifies how many DOFs are associated with an entity’s interior.

### Local to Global Index Mapping
DOF handlers use index mapping to associate local shape functions on elements to global DOFs. The mapping follows a convention to ensure consistent assembly of the global system. The basic function is represented as:

```cpp
locglobmap(K, i) = j
```

Where:

- K refers to a mesh entity (such as a cell).
- i is the local index of the shape function.
- j is the corresponding global index.
- For example, for linear Lagrangian elements on triangular meshes, the following local-to-global map might apply:

```cpp
locglobmap(K, 0) = 2
locglobmap(K, 1) = 7
locglobmap(K, 2) = 9
```

Global shape functions (DOFs) in LehrFEM++ follow a numbering convention based on the geometric entities they are associated with:

Points (Vertices)
1. Edges
2. Cells (Triangles or Quadrilaterals)
3. DOFs associated with lower-dimensional entities (e.g., points and edges) are numbered first, followed by higher-dimensional ones (cells). Within each category, entities are numbered according to their intrinsic indexing as returned by the Index() function.

### Example: UniformFEDofHandler

A commonly used DOF handler in LehrFEM++ is the lf::assemble::UniformFEDofHandler. It sets up local-to-global index mappings where the number of DOFs is constant for all entities of the same type:

```cpp
lf::assemble::UniformFEDofHandler dof_handler(
  mesh_p, {{lf::base::RefEl::kPoint(), 1},
           {lf::base::RefEl::kSegment(), 2},
           {lf::base::RefEl::kTria(), 1},
           {lf::base::RefEl::kQuad(), 4}});
```
This means:

- Each node carries 1 DOF.
- Each edge carries 2 DOFs.
- Each triangle has 1 DOF.
- Each quadrilateral has 4 DOFs.

### Special Conventions
In LehrFEM++, DOFs associated with lower-dimensional entities are numbered first. For example, DOFs on points (nodes) are numbered before those on edges. This is simply a convention to ensure consistent assembly and output.

### Important Methods
Here are some common functions and their usage in DofHandler:

- DOF Information Extraction: The GlobalDofIndices() method can be used to extract global DOF indices covering any geometric entity. Similarly, NumLocalDofs() provides the count of DOFs on each entity.
- Interior DOFs: The InteriorGlobalDofIndices() function returns the indices of global shape functions associated with the interior of a given entity. This is particularly useful when dealing with higher-order elements where some DOFs are attached to edges or faces rather than vertices.

### Example Code
An example of using a DofHandler to print information about local-to-global mappings in LehrFEM++:

```cpp
void printDofInfo(const lf::assemble::DofHandler &dofh) {
  auto mesh = dofh.Mesh();
  size_type N_dofs = dofh.NumDofs();
  
  std::cout << "DofHandler (" << N_dofs << " dofs):\n";
  
  for (lf::base::dim_t codim = 0; codim <= mesh->DimMesh(); ++codim) {
    for (const lf::mesh::Entity *e : mesh->Entities(codim)) {
      size_type no_dofs = dofh.NumLocalDofs(*e);
      auto dof_array = dofh.GlobalDofIndices(*e);
      
      std::cout << *e << " : " << no_dofs << " dofs = [";
      for (int loc_dof_idx = 0; loc_dof_idx < no_dofs; ++loc_dof_idx) {
        std::cout << dof_array[loc_dof_idx] << " ";
      }
      std::cout << "]\n";
    }
  }
}
```
