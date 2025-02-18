# Quick Reference - Boundary Conditions {#quick_reference_bc}

[TOC]

> [!caution]
> The contents of this page are discussed in @lref_link{sec:essbdc}. Please read this section before using the quick reference.

## Overview

LehrFEM++ provides a number of convenient functions for working with essential boundary conditions.

To fix degrees of freedom (DoFs) on the boundary, we can use the functions `lf::assemble::FixFlaggedSolutionComponents` or `lf::assemble::FixFlaggedSolutionCompAlt` (see also [Quick Reference - Assembly](@ref quick_reference_assembly)). Both functions take a function as an argument that returns a pair of a boolean and a double: `std::pair<bool, double> selector(unsigned int dof_idx)`. The selector function returns whether the DoF is to be fixed and, if so, the value it should be fixed to.

We can obtain the boundary flags as a [MeshDataSet](@ref mesh_data_sets) using:

```cpp
auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2);
```

The second argument is the co-dimension of the entities we are interested in. To fix all DoFs on the boundary (including those associated with edges and cells), omit the second argument.

## One function on the whole boundary {#one_function}

Let \f$g\f$ be a function defining the values on the boundary. To get the function values of \f$g\f$ at the boundary nodes, wrap it in a MeshFunction:

```cpp
auto mf_g = lf::mesh::utils::MeshFunctionGlobal(g);

std::vector<std::pair<bool, scalar_t>> boundary_val = 
    lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g);
```

This returns a `std::vector<std::pair<bool, scalar_t>>`. To create our selector:

```cpp
auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
    return boundary_val[dof_idx];
};
```

If \f$g\f$ has different definitions on different parts of the boundary (e.g., \f$g = 0\f$ on \f$\gamma_0\f$ and \f$g = 1\f$ on
\f$\gamma_1\f$), try to express \f$g\f$ as a lambda function with an if-else statement.

## Constant value on the boundary {#constant_value}

The selector becomes:

```cpp
auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
    if (bd_flags[dof_idx]) {
        return std::make_pair(true, boundary_value);
    } else {
        return std::make_pair(false, 0.0); // value irrelevant
    }
};
```

## BC only on part of the boundary {#part_of_boundary}

We need to create our own `bd_flags`. Initialize a `MeshDataSet` with default value `false`:

```cpp
lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(mesh_p, false);
```

Then loop over nodes and edges separately and set `bd_flags` to true for the nodes/edges where the boundary condition should be applied.

```cpp
for (const auto& edge : fe_space->Mesh()->Entities(1)) {
    if (...) bd_flags(*edge) = true;
}

for (const auto& node : fe_space->Mesh()->Entities(2)) {
    // 
    ...
}
```

