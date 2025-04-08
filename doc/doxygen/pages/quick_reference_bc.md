# Quick Reference - Boundary Conditions {#quick_reference_bc}

[TOC]

@gh_edit

> [!caution]
> The contents of this page are discussed in @lref_link{sec:essbdc}. Please read this section before using the quick reference.

## Overview

LehrFEM++ provides a number of convenient functions for working with essential boundary conditions.

To fix degrees of freedom (DoFs) on the boundary, we can use the functions `lf::assemble::FixFlaggedSolutionComponents` or `lf::assemble::FixFlaggedSolutionCompAlt` (see also [Quick Reference - Assembly](@ref quick_reference_assembly)). Both functions take a function as an argument that returns a pair of a boolean and a double: `std::pair<bool, double> selector(unsigned int dof_idx)`. The selector function returns whether the DoF is to be fixed and, if so, the value it should be fixed to.

We can obtain the boundary flags as a [MeshDataSet](@ref mesh_data_sets) using:

@snippet{trimleft} quick_reference/qr_bc_snippets.cpp bc Overview

The second argument is the co-dimension of the entities we are interested in. To fix all DoFs on the boundary (including those associated with edges and cells), omit the second argument.

## One function on the whole boundary {#one_function}

Let \f$g\f$ be a function defining the values on the boundary. To get the function values of \f$g\f$ at the boundary nodes, wrap it in a MeshFunction:

@snippet{trimleft} quick_reference/qr_bc_snippets.cpp bc One function on the whole boundary

This returns a `std::vector<std::pair<bool, scalar_t>>`. To create our selector:

@snippet{trimleft} quick_reference/qr_bc_snippets.cpp bc One function on the whole boundary 0

If \f$g\f$ has different definitions on different parts of the boundary (e.g., \f$g = 0\f$ on \f$\gamma_0\f$ and \f$g = 1\f$ on
\f$\gamma_1\f$), try to express \f$g\f$ as a lambda function with an if-else statement.

## Constant value on the boundary {#constant_value}

The selector becomes:

@snippet{trimleft} quick_reference/qr_bc_snippets.cpp bc Constant value on the boundary

## BC only on part of the boundary {#part_of_boundary}

We need to create our own `bd_flags`. Initialize a `MeshDataSet` with default value `false`:

@snippet{trimleft} quick_reference/qr_bc_snippets.cpp bc BC only on part of the boundary

Then loop over nodes and edges separately and set `bd_flags` to true for the nodes/edges where the boundary condition should be applied.

@snippet{trimleft} quick_reference/qr_bc_snippets.cpp bc BC only on part of the boundary 0


## Example (Essential boundary conditions) {#example}

@snippet{trimleft} fe_tools.cc InitEssentialConditionFromFunction
