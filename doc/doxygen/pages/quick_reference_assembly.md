# Quick Reference - Assembly {#quick_reference_assembly}

[TOC]

@gh_edit

> [!caution]
> The contents of this page are discussed in @lref_link{sec:assembly}. Please read this section before using the quick reference.

## Overview

Assembly in the context of LehrFEM++ refers to computing entries of the stiffness matrix and the right-hand side vector (load vector) for finite elements. @lref_link{sec:assembly} provides more detailed information.

A key part of efficient assembly is cell-local computations, where the indexing is handled by [DOFHandlers](@ref quick_reference_dofh). To solve a PDE using the Galerkin method, we _assemble_ the matrix \f$\textbf{A}\f$ and the right-hand side vector \f$\vec{b}\f$ to obtain a linear system \f$\textbf{A}\vec{x} = \vec{b}\f$.

LehrFEM++ provides a number of helper functions for assembly, which are detailed in @lref_link{sss:assalg}.

## Stiffness/Galerkin Matrix

LehrFEM++ provides the `lf::assemble::AssembleMatrixLocally` function for assembling the stiffness matrix \f$\textbf{A}\f$. The following example shows how to assemble the stiffness matrix for a Poisson bilinear form

\f[
a(u, v) = \int_{\Omega} \textbf{grad} u \cdot \textbf{grad} v \, dx
\f]

using the [entity matrix provider](@ref quick_reference_emp) `lf::uscalfe::LinearFELaplaceElementMatrix`.

@snippet snippets/assembler.cc matrix_usage

## Right-Hand Side (RHS)/Load Vector

LehrFEM++ also provides the `lf::assemble::AssembleVectorLocally` function for assembling the right-hand side vector \f$\vec{b}\f$. The following example shows how to assemble a right-hand side vector with a custom [entity vector provider](@ref quick_reference_emp), where \f$b\f$ is assembled from volume contributions from local cells.

@snippet snippets/assembler.cc vector_usage