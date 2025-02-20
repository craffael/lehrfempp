# Quick Reference - Quadrature {#quick_reference_quad}

[TOC]

@gh_edit

> [!caution]
> The contents of this page is discussed in @lref_link{par:lfquad}. Please read this section before using the quick reference.

# Overview

LehrFEM++ provides a number of quadrature rules for the numerical integration of functions over reference elements. General entities can be mapped to the reference element using the methods described in the [Geometry](@ref quick_reference_geometry) quick reference.

# Working with Quadrature Rules
The `quad` namespace has some useful functions for instantiating quadrature rules for a given reference element. The following code snippet shows how to get a quadrature rule for a triangle of degree >= 3 using `lf::quad::make_QuadRule`.

@snippet{trimleft} quick_reference/qr_quad_snippets.cpp quad Working with Quadrature Rules

It is also possible to get a quadrature rule which only uses the nodes of the reference element to approximate the integral using `lf::quad::make_QuadRuleNodal`.

@snippet{trimleft} quick_reference/qr_quad_snippets.cpp quad Working with Quadrature Rules 0

Another possibility is to define a quadrature rule by specifying the nodes and weights manually. This can be done using the `lf::quad::QuadRule` constructor.

@snippet{trimleft} quick_reference/qr_quad_snippets.cpp quad Working with Quadrature Rules 1

Lastly LehrFEM++ also provides a few more specialized quadrature rules. A list can be found on the [Quadrature Rules Namespace](@ref lf::quad) page.

# Using a Quadrature Rule

@snippet{trimleft} quick_reference/qr_quad_snippets.cpp quad Using a Quadrature Rule

