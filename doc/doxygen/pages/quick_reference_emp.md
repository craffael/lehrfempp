# Quick Reference - Enitiy Matrix Providers {#quick_reference_emp}

[TOC]

@gh_edit

> [!caution]
> The contents of this page is discussed in @lref_link{rem:cppconcepts}. Please read this section before using the quick reference.

## Overview

LehrFEM++ includes built-in `lf::assemble::EntityMatrixProvider` (EMPs) and `lf::assemble::EntityVectorProvider` (EVPs) for common PDEs. Custom EMPs/EVPs can also be defined. These providers compute local element matrices/vectors for finite element methods and are typically passed to `lf::assemble::AssembleMatrixLocally` or `lf::assemble::AssembleVectorLocally` (see [Assembly](@ref quick_reference_assembly)).

## Custom Entity Vector Providers {#custom_entity_vector_providers}

A custom EMP must implement the `lf::assemble::EntityMatrixProvider` concept. Minimal definitions:

@snippet snippets/assembler.cc lflinfeelmat

@snippet snippets/assembler.cc lflinfeelvec

## Built-in Entity Matrix Providers {#built_in_entity_matrix_providers}

<!-- Editor Note: Turn off word wrap to display table properly in Editor -->

Below is a summary of key built-in EMPs:

| Equation                                                                                                                                                                                              | Name (follow link for details)                                                  | Description                                           |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- | ----------------------------------------------------- |
| \f$ \int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u \cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x} \f$                                                                                    | [fe::DiffusionElementMatrixProvider](#lffediffusionelementmatrixprovider)       | stiffness matrix for the diffusion equation.          |
| \f$ \int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x} \f$                                                                                                                        | [fe::MassElementMatrixProvider](#lfemasselementmatrixprovider)                  | mass matrix for the reaction-diffusion equation.      |
| \f$ \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x) \f$                                                                                                                                     | [fe::MassEdgeMatrixProvider](#lfemassedgematrixprovider)                        | mass matrix for the reaction-diffusion equation.      |
| \f$ \int\limits_{K}\mathbf{grad}\,u\cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x} \f$                                                                                                                    | [uscalfe::LinearFELaplaceElementMatrix](#lfuscalfelinearfelaplaceelementmatrix) | stiffness matrix for the Laplace equation.            |
| \f$ \int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u\cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x} \f$ \n \f$ + \int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x} \f$ | [uscalfe::ReactionDiffusionElementMatrixProvider](#lfuscalferdemp)              | stiffness matrix for the reaction-diffusion equation. |
| \f$ \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x) \f$                                                                                                                                     | [uscale::MassEdgeMatrixProvider](#lfuscalfemassedgematrixprovider)              | mass matrix for the reaction-diffusion equation.      |

### lf::fe::DiffusionElementMatrixProvider {#lffediffusionelementmatrixprovider}

lf::fe::DiffusionElementMatrixProvider

\f[
    (u,v) \mapsto\int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u
          \cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x}
 \;
\f]

with diffusion coefficient \f$\mathbf{\alpha}\f$, see also Example @lref_link{ex:rdemp}

### lf::fe::MassElementMatrixProvider {#lfemasselementmatrixprovider}

lf::fe::MassElementMatrixProvider

The element matrix corresponds to the (local) bilinear form
\f[
    (u,v)
 \mapsto\int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x}
 \;,
\f]

with reaction coefficient \f$\gamma\f$, see also @lref_link{ex:rdemp}

### lf::fe::MassEdgeMatrixProvider {#lfemassedgematrixprovider}

lf::fe::MassEdgeMatrixProvider

The edge matrix corresponds to the (local) bilinear form

\f[
    (u,v) \mapsto \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x)\;,
\f]

where @f$e@f$ is an edge of the mesh, and \f$\gamma\f$ a scalar-valued coefficient function. 

### lf::uscalfe::LinearFELaplaceElementMatrix {#lfuscalfelinearfelaplaceelementmatrix}

lf::uscalfe::LinearFELaplaceElementMatrix

The element matrix corresponds to the (local) bilinear form

\f[
    (u,v)
 \mapsto\int\limits_{K}\mathbf{grad}\,u\cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x}
 \;,
\f]


### lf::uscalfe::ReactionDiffusionElementMatrixProvider {#lfuscalferdemp}

lf::uscalfe::ReactionDiffusionElementMatrixProvider


\f[
    (u,v) \mapsto\int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u\cdot\mathbf{grad}\,v +
    \gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x}
 \;,
\f]

with diffusion coefficient \f$\mathbf{\alpha}\f$ and reaction coefficient \f$\gamma\f$. See also @lref_link{ex:rdemp}.

### lf::uscalfe::MassEdgeMatrixProvider {#lfuscalfemassedgematrixprovider}

lf::uscalfe::MassEdgeMatrixProvider

\f[
    (u,v) \mapsto \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x)\;,
\f]

where @f$e@f$ is an edge of the mesh, and \f$\gamma\f$ a scalar-valued coefficient function. See @lref_link{ex:lfemp} for a similar example on assembly of boundary contributions.

## Built-in Entity Vector Providers

Entity Vector Providers (EVPs) in LehrFEM++ are used to compute local element vectors for finite element methods. Below is a table summarizing the key built-in EVPs in LehrFEM++.

| Equation                                                                  | Name (follow link for details)                                                        | Description                                                            |
| ------------------------------------------------------------------------- | ------------------------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| \f$ \int\limits_{K}f(\mathbf{x})v\,\mathrm{d}\mathbf{x} \f$               | [fe::ScalarLoadElementVectorProvider](#lffeScalarLoadElementVectorProvider)           | load vector for the scalar load on an element.                         |
| \f$ \int\limits_{e}g(x)\overline{v(x)}\,\mathrm{d}S(x) \f$                | [fe::ScalarLoadEdgeVectorProvider](#lffeScalarLoadEdgeVectorProvider)                 | load vector for the scalar load on an edge \f$e\f$.                    |
| \f$ \int\limits_{K}f(\mathbf{x})\overline{v(x)}\,\mathrm{d}\mathbf{x} \f$ | [uscalfe::ScalarLoadElementVectorProvider](#lfuscalfeScalarLoadElementVectorProvider) | load vector for scalar finite elements; **volume contributions only**. |
| \f$ \int\limits_{e}g(x)\overline{v(x)}\,\mathrm{d}S(x) \f$                | [uscalfe::ScalarLoadEdgeVectorProvider](#lfuscalfeScalarLoadEdgeVectorProvider)       | load vector for scalar finite elements                                 |
| \f$ \int\limits_{K}f(\mathbf{x})\overline{v(x)}\,\mathrm{d}\mathbf{x} \f$ | [uscalfe::LinearFELocalLoadVector](#lfuscalfeLinearFELocalLoadVector)                 | load vector for scalar **linear** finite elements.                     |

### lf::fe::ScalarLoadElementVectorProvider {#lffeScalarLoadElementVectorProvider}

lf::fe::ScalarLoadElementVectorProvider computes the local element vector corresponding to the (local) linear form

\f[
    v \mapsto\int\limits_{K}f(\mathbf{x})v\,\mathrm{d}\mathbf{x}
 \;,
\f]

where \f$f\f$ is a scalar-valued function.

### lf::fe::ScalarLoadEdgeVectorProvider {#lffeScalarLoadEdgeVectorProvider}

lf::fe::ScalarLoadEdgeVectorProvider computes the local edge contributions corresponding to the (local) linear form

\f[
    v \mapsto\int\limits_{e}g(x)\overline{v(x)}\,\mathrm{d}S(x)
 \;,
\f]

where \f$e\f$ is an edge of the mesh, and \f$g\f$ a scalar-valued function.

### lf::uscalfe::ScalarLoadElementVectorProvider {#lfuscalfeScalarLoadElementVectorProvider}

lf::uscalfe::ScalarLoadElementVectorProvider computes the local element vector corresponding to the (local) linear form

\f[
    v \mapsto\int\limits_{K}f(\mathbf{x})\overline{v(x)}\,\mathrm{d}\mathbf{x}
 \;,
\f]

where \f$f\f$ is a scalar-valued function. **Note**: This provider only computes the volume contributions on \f$K\f$.

### lf::uscalfe::ScalarLoadEdgeVectorProvider {#lfuscalfeScalarLoadEdgeVectorProvider}

lf::uscalfe::ScalarLoadEdgeVectorProvider computes the local edge contributions corresponding to the (local) linear form

\f[
    v \mapsto\int\limits_{e}g(x)\overline{v(x)}\,\mathrm{d}S(x)
 \;,
\f]

where \f$e\f$ is an edge of the mesh, and \f$g\f$ a scalar-valued function.

### lf::uscalfe::LinearFELocalLoadVector {#lfuscalfeLinearFELocalLoadVector}

lf::uscalfe::LinearFELocalLoadVector computes the local element vector for **linear** FE corresponding to the (local) linear form 

\f[
    v \mapsto\int\limits_{K}f(\mathbf{x})\overline{v(x)}\,\mathrm{d}\mathbf{x}
 \;,
\f]

where \f$f\f$ is a scalar-valued function using edge midpoint quadrature.


<!-- ## Example

@snippet examples/lecturedemos/lecturedemoassemble.cc  -->
