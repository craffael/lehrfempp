# Quick Reference - Boundary Conditions {#quick_reference_bc}

[TOC]

LehrFEM provides a number of convenient functions to work with essential boundary conditions. More information can be found in the [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf) @lref{sec:essbdc}.

To fix degrees of freedom on the boundary we can use the functions lf::assemble::FixFlaggedSolutionComponents or lf::assemble:: FixFlaggedSolutionCompAlt (see also [QR Assembly](./quick_reference_assembly.md)). Both functions takes a function as an argument that returns a pair of a boolean and a double: `std::pair<bool, double> selector(unsigned int dof_idx)`. The selector function returns whether the dof is to be fixed and if so, the
value it should be fixed to.

We get the boundary flags as a [MeshDataSet](./quick_reference_mesh.html#mesh_data_sets) using:

```cpp
auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2);
```

The second argument is the co-dimension of entities we are interested in. To fix all [dofs](./quick_reference_dofs.html) on the boundary (including those associated with edges and cells), give no second argument. 

## One function on the whole boundary {#one_function}

Let \f$g\f$ be a function defining the values on the boundary. To get the function values of \f$g\f$ at the boundary nodes, wrap it in a MeshFunction:

```cpp
auto mf_g = lf::mesh::utils::MeshFunctionGlobal(g);
auto boundary_val = lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g);
```

This gives a std::vector<std::pair<bool, scalar_t>>. To get our selector:

```cpp
auto selector = [&](unsigned int dof_idx) -> std::pair<bool, double> {
    return boundary_val[dof_idx];
};
```

If \f$g\f$ has different definitions on different parts of the boundary (e.g., \f$g = 0\f$ on \f$\gamma_0\f$ and \f$g = 1\f$ on
\f$\gamma_1\f$), try to express \f$g\f$ as a lambda function with an if-else statement.

## Constant value on the boundary {#constant_value}

The selector becomes
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

Wee need to create our own bd_flags. Initialize a MeshDataSet with default value `false`:

```cpp
lf::mesh::utils::AllCodimMeshDataSet<bool> bd_flags(mesh_p, false);
```

Then loop over nodes and edges separately and set `bd_flags` to true for the nodes/edges where
the BC should be applied.

```cpp
for (const auto& edge : fe_space->Mesh()->Entities(1)) {
    if (...) bd_flags(*edge) = true;
}

for (const auto& node : fe_space->Mesh()->Entities(2)) {
    // 
    ...
}
```

