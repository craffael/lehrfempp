# Quick Reference - Enitiy Matrix Providers {#quick_reference_emp}

[TOC]

> [!caution]
> The contents of this page is discussed in @lref_link{rem:cppconcepts}. Please read before using quick reference.

## Overview

LehrFEM++ provides some built-in lf::assemble::EntityMatrixProvider (EMPs) and entity lf::assemble::EntityVectorProvider (EVPs) for common PDEs. Users can also define custom EMPs and EVPs. EMPs/EVPs are used to compute local element matrices for finite element methods and are usually passed to the `lf::assemble::AssembleMatrixLocally` / `lf::assemble::AssembleVectorLocally` functions (see also [Quick Reference - Assembly](@ref quick_reference_assembly)).

## Custom Entity Vector Providers {#custom_entity_vector_providers}

A custom entity matrix provider has to implement the `lf::assemble::EntityMatrixProvider` concept. The following snippets offer minimal definitions for EMPs and EVPs:

@snippet snippets/assembler.cc lflinfeelmat

@snippet snippets/assembler.cc lflinfeelvec

## Built-in Entity Matrix Providers {#built_in_entity_matrix_providers}

<!-- TODO (barmstron): Write down quick specification of each EL/VE MP -->

<!-- Editor Note: Turn off word wrap to display table properly in Editor -->

Entity Matrix Providers (EMPs) in LehrFEM++ are used to compute local element matrices for finite element methods. Below is a table summarizing the key built-in EMPs in LehrFEM++.

| Equation                                                                                                                                                                        | Name (follow link for details)                                                  | Description                                                              |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------- | ------------------------------------------------------------------------ |
| \f$ \int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u \cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x} \f$                                                              | [fe::DiffusionElementMatrixProvider](#lffediffusionelementmatrixprovider)       | Computes the local stiffness matrix for the diffusion equation.          |
| \f$ \int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x} \f$                                                                                                  | [fe::MassElementMatrixProvider](#lfemasselementmatrixprovider)                  | Computes the local mass matrix for the reaction-diffusion equation.      |
| \f$ \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x) \f$                                                                                                               | [fe::MassEdgeMatrixProvider](#lfemassedgematrixprovider)                        | Computes the local mass matrix for the reaction-diffusion equation.      |
| \f$ \int\limits_{K}\mathbf{grad}\,u\cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x} \f$                                                                                              | [uscalfe::LinearFELaplaceElementMatrix](#lfuscalfelinearfelaplaceelementmatrix) | Computes the local stiffness matrix for the Laplace equation.            |
| \f$ \int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u\cdot\mathbf{grad}\,v + \f$ \n \f$ \int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x} \f$ | [uscalfe::ReactionDiffusionElementMatrixProvider](#lfuscalferdemp)              | Computes the local stiffness matrix for the reaction-diffusion equation. |
| \f$ \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x) \f$                                                                                                               | [uscale::MassEdgeMatrixProvider](#lfuscalfemassedgematrixprovider)              | Computes the local mass matrix for the reaction-diffusion equation.      |

### lf::fe::DiffusionElementMatrixProvider {#lffediffusionelementmatrixprovider}

lf::fe::DiffusionElementMatrixProvider

$$
    (u,v) \mapsto\int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u
          \cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x}
 \;
$$

with diffusion coefficient \f$\mathbf{\alpha}\f$, see also Example [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf) @lref{ex:rdemp}

### lf::fe::MassElementMatrixProvider {#lfemasselementmatrixprovider}

lf::fe::MassElementMatrixProvider

The element matrix corresponds to the (local) bilinear form
$$
    (u,v)
 \mapsto\int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x}
 \;,
$$

with reaction coefficient \f$\gamma\f$, see also [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf) @lref{ex:rdemp}

### lf::fe::MassEdgeMatrixProvider {#lfemassedgematrixprovider}

lf::fe::MassEdgeMatrixProvider

The edge matrix corresponds to the (local) bilinear form

$$
    (u,v) \mapsto \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x)\;,
 $$

where @f$e@f$ is an edge of the mesh, and \f$\gamma\f$ a scalar-valued coefficient function.

### lf::uscalfe::LinearFELaplaceElementMatrix {#lfuscalfelinearfelaplaceelementmatrix}

lf::uscalfe::LinearFELaplaceElementMatrix

The element matrix corresponds to the (local) bilinear form

$$
    (u,v)
 \mapsto\int\limits_{K}\mathbf{grad}\,u\cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x}
 \;,
$$


### lf::uscalfe::ReactionDiffusionElementMatrixProvider {#lfuscalferdemp}

lf::uscalfe::ReactionDiffusionElementMatrixProvider


$$
    (u,v) \mapsto\int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u\cdot\mathbf{grad}\,v +
    \gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x}
 \;,
$$

### lf::uscalfe::MassEdgeMatrixProvider {#lfuscalfemassedgematrixprovider}

lf::uscalfe::MassEdgeMatrixProvider

$$
    (u,v) \mapsto \int\limits_e\gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x)\;,
$$

where @f$e@f$ is an edge of the mesh, and \f$\gamma\f$ a scalar-valued coefficient function.

## Built-in Entity Vector Providers

Entity Vector Providers (EVPs) in LehrFEM++ are used to compute local element vectors for finite element methods. Below is a table summarizing the key built-in EVPs in LehrFEM++.

| Equation                                                    | Name (follow link for details)                                             | Description                                                       |
| ----------------------------------------------------------- | -------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| \f$ \int\limits_{K}f(\mathbf{x})v\,\mathrm{d}\mathbf{x} \f$ | [fe::ScalarLoadElementVectorProvider](#lfescalarloadelementvectorprovider) | Computes the local load vector for the scalar load on an element. |

### lf::fe::ScalarLoadElementVectorProvider {#lfescalarloadelementvectorprovider}
### lf::fe::ScalarLoadEdgeVectorProvider {#lfescalarloadedgevectorprovider}
### lf::uscalfe::LinearFELocalLoadVector {#lfuscalfelelementvectorprovider}
### lf::uscalfe::ScalarLoadElementVectorProvider {#lfuscalfeelementvectorprovider}
### lf::uscalfe::ScalarLoadEdgeVectorProvider {#lfuscalfeedgevectorprovider}
