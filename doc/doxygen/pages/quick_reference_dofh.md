# Quick Reference - Degrees of Freedom (DOF) Handlers and Indexing {#quick_reference_dofh}

[TOC]

@gh_edit

> [!caution]
> The contents of this page is discussed in @lref_link{sec:parmFE}. Please read this section before using the quick reference.
## Overview

Local computations are essential for efficient finite element methods (FEM). LehrFEM++ uses a local indexing scheme to map local shape functions to global degrees of freedom (DOFs). In **LehrFEM++**, this is done through objects of the type `lf::assemble::DofHandler`. This quick reference provides an overview of DOF handlers and indexing conventions in **LehrFEM++**. The [assembly quick reference](@ref quick_reference_assembly) page provides more information on (local) assembly using DOF handlers.
<!-- In finite element methods (FEM), handling degrees of freedom (DOFs) is critical for constructing systems of equations. DOF handlers are responsible for assigning local basis functions to global basis functions, which allows for the assembly of global matrices and vectors from element contributions. -->

## Local to Global Index Mapping

The local shape functions of every cell (co-dimension-0 entity) are arranged according to the increasing dimension of the geometric entities they are associated with:

1. Points (Vertices)
2. Edges
3. Cells (Triangles or Quadrilaterals)

DOFs associated with lower-dimensional entities (e.g., points and edges) are numbered first, followed by higher-dimensional ones (cells). Within each category, entities are numbered according to their intrinsic indexing as returned by the `Index()` function. This is also illustrated in @lref_link{quadnodes2} and @lref_link{lnlfe2}.

## DOF Handlers in LehrFEM++

LehrFEM++ provides two implementations of DOF handlers:

1. `lf::assemble::DynamicFEDofHandler`: Allows for a variable number of local DOFs for each entity (e.g. for hp-FEM, see @lref_link{par:dofhinit})
2. `lf::assemble::UniformFEDofHandler`: Assigns the same number of DOFs to each entity of a given type (uniform FE space).


### DynamicFEDofHandler

A `DynamicFEDofHandler` can be initialized with a function that returns the number of DOFs associated with a given entity.

```cpp
const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

auto func = [](const lf::mesh::Entity& e) {
  // return number of dofs associated with e
  return 1;
};

lf::assemble::DynamicFEDofHandler dofh(mesh_p, func);
```

### UniformFEDofHandler

```cpp
// Generate a simple test mesh
const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
    lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

// Initialize a DofHandler for p=2 Lagrange FE
lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                        {{lf::base::RefEl::kPoint(), 1},
                                        {lf::base::RefEl::kSegment(), 1},
                                        {lf::base::RefEl::kTria(), 0},
                                        {lf::base::RefEl::kQuad(), 1}});
```

For example, in a second-order Lagrange FE space (as shown in the code above):

- Each point carries 1 DOF.
- Each edge carries 1 DOF.
- Each triangle has 0 DOFs.
- Each quadrilateral has 1 DOF.

## Global Indexing convention

LehrFEM++ follows the convention of numbering **global** DOFs according to the following rule:

1. DOFs associated with lower-dimensional entities are numbered first:
   \f[
    \text{POINT} \rightarrow \text{SEGMENT} \rightarrow \text{TRIA/QUAD}
   \f]

2. The indices of DOFs belonging to entities of the same co-dimension increase with increasing entity indices as returned by the `Index()` function.

## Key Methods
The key methods provided by the `DofHandler` interface are:

- `NumDofs()`: Returns the total number of global DOFs. See also @lref_link{par:betldofmap}.

  ```cpp
  unsigned num_dofs = dofh.NumDofs();
  ```

- `NumLocalDofs(const lf::mesh::Entity &)`: Returns the number of DOFs associated with a particular geometric entity.  See also @lref_link{par:betldofmap}.
  
  ```cpp
  const lf::mesh::Entity* e = mesh_p->EntityByIndex(0, 0);

  unsigned num_local_dofs = dofh.NumLocalDofs(*e);
  ```
- `GlobalDofIndices(const lf::mesh::Entity &)`: Returns the global indices of DOFs associated with a given entity.  See also @lref_link{par:betldofmap}.

  ```cpp
  const lf::mesh::Entity* e = mesh_p->EntityByIndex(0, 0);

  std::span<const lf::assemble::gdof_idx_t> idx = dof2.GlobalDofIndices(*e);

  for (auto i : idx) {
    std::cout << i << " ";
  } // -> prints the global indices of DOFs associated with the entity

  ```

- `NumInteriorDofs(const lf::mesh::Entity &)`: Specifies how many DOFs are associated with an entity’s interior.  See also @lref_link{par:betldofmap}.

  ```cpp
  const lf::mesh::Entity* e = mesh_p->EntityByIndex(0, 0);

  unsigned num_int_dofs = dofh.NumInteriorDofs(*e);
  ```

## Example Code
An example of using a DofHandler to print information about local-to-global mappings in LehrFEM++:

@snippet snippets/dofhandleruse.cc pdi
