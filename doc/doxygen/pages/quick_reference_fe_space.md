# Quick Reference - Finite Element Spaces {#quick_reference_fe_space}

[TOC]

@gh_edit

> [!caution]
> The contents of this page are discussed in @lref_link{par:fespace} and @lref_link{sec:LagrFEM}. Please read before using quick reference.

## Overview

LehrFEM++ provides an interface for (scalar) finite element spaces. Scalar means that the approximation space is always scalar-valued, and the shape functions are scalar-valued.

## (General) Scalar Finite Element Space (lf::fe) {#finite_element_space}

The `lf::fe` namespace provides a general interface for scalar finite element spaces. Implementation and usage of scalar FE spaces are covered in more detail in @lref_link{par:fespace}.

### Hierarchic Scalar Finite Element Space

A special case of the general scalar finite element space is a **Hierarchic Scalar Finite Element Space**, implemented by `lf::fe::HierarchicScalarFESpace`. Hierarchic FE spaces assign a polynomial degree for the shape functions to each mesh entity. They can easily be constructed from a mesh and a function mapping mesh entities to polynomial degrees.

@snippet{trimleft} quick_reference/qr_fe_space_snippets.cpp fe_space Hierarchic Scalar Finite Element Space

## Uniform Finite Element Space (lf::uscalfe) {#uniform_finite_element_spaces}

The `lf::uscalfe` (<strong>U</strong>niform <strong>Scal</strong>ar <strong>F</strong>inite <strong>E</strong>lements) namespace is a specialization of `lf::fe` for uniform scalar finite element spaces. Uniform in this context means that the approximation space has a uniform order of approximation over the whole mesh. In other words, the shape functions of a given approximation space depend only on the underlying reference element of a mesh entity.

### Lagrangian Finite Element Spaces {#fe_space_lagrange}

A prominent example of a uniform finite element space is the n-th order Lagrangian finite element spaces. @lref_link{sec:LagrFEM} discusses these spaces and the mathematical background in detail.

LehrFEM++ provides convenience classes for constructing order 1, 2, and 3 Lagrangian finite element spaces:

- `lf::uscalfe::FeSpaceLagrangeO1`: (Bi)Linear Lagrangian Finite Element space.
- `lf::uscalfe::FeSpaceLagrangeO2`: Quadratic Lagrangian Finite Element space.
- `lf::uscalfe::FeSpaceLagrangeO3`: Cubic Lagrangian Finite Element space.

@snippet{trimleft} quick_reference/qr_fe_space_snippets.cpp fe_space Lagrangian Finite Element Spaces
