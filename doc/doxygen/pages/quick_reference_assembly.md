# Quick Reference - Assembly {#quick_reference_assembly}

[TOC]

> [!caution]
> The contents of this page is discussed in @lref_link{sec:assembly}. Please read before using quick reference.

## Overview

Assembly in the context of LehrFEM++ refers to computing entries of stiffness matrix/right hand side vector (load vector) for finite elements. @lref_link{sec:assembly} goes into detail.

A key part of efficient assembly are cell-local computations for where the indexing is handled by [DOFHandlers](@ref quick_reference_dofh). To solve a PDE using the Galerkin method we _assemble_ the matrix \f$\textbf{A}\f$ and the right hand side vector \f$\vec{b}\f$ to obtain a linear system \f$\textbf{A}\vec{x} = \vec{b}\f$.

LehrFEM++ provides a number of helper functions for assembly detailed in @lref_link{sss:assalg}.

## Stiffness/Galerkin matrix

LehrFEM++ provides the `lf::assemble::AssembleMatrixLocally` function for assembling the stiffness matrix \f$\textbf{A}\f$. The following example shows how to assemble the stiffness Matrix for a poisson bi-linear form

\f[
a(u, v) = \int_{\Omega} \textbf{grad} u \cdot \textbf{grad} v \, dx
\f]

through the [entity matrix provider](@ref quick_reference_emp) `lf::uscalfe::LinearFELaplaceElementMatrix`.

@snippet snippets/assembler.cc matrix_usage

## Right hand side (RHS)/load vector

LehrFEM++ also provides the `lf::assemble::AssembleVectorLocally` function for assembling the right hand side vector \f$\vec{b}\f$. The following example shows how to assemble a right hand side vector with a custom [entity vector provider](@ref quick_reference_emp) where \f$b\f$ is assembled from volume contributions from local cells.

@snippet snippets/assembler.cc vector_usage

