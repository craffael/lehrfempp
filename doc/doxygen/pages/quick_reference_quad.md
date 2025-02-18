# Quick Reference - Quadrature {#quick_reference_quad}

[TOC]

> [!caution]
> The contents of this page is discussed in @lref_link{par:lfquad}. Please read this section before using the quick reference.

# Overview

LehrFEM++ provides a number of quadrature rules for the numerical integration of functions over reference elements. General entities can be mapped to the reference element using the methods described in the [Geometry](@ref quick_reference_geometry) quick reference.

# Working with Quadrature Rules
The `quad` namespace has some useful functions for instantiating quadrature rules for a given reference element. The following code snippet shows how to get a quadrature rule for a triangle of degree >= 3 using `lf::quad::make_QuadRule`.

```cpp
auto ref_tria = lf::base::RefEl::kTria();

lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_tria, 3);
```

It is also possible to get a quadrature rule which only uses the nodes of the reference element to approximate the integral using `lf::quad::make_QuadRuleNodal`.

```cpp
lf::quad::QuadRule qr = lf::quad::make_QuadRuleNodal(ref_tria);
```

Another possibility is to define a quadrature rule by specifying the nodes and weights manually. This can be done using the `lf::quad::QuadRule` constructor.

```cpp
Eigen::MatrixXd points(2, 3);
points << 0.166667, 0.666667, 0.166667,
          0.166667, 0.166667, 0.666667;

Eigen::Vector3d weights;
weights << 0.166667, 0.166667, 0.166667;

auto qr = lf::quad::QuadRule(ref_tria, points, weights, 2);
```

Lastly LehrFEM++ also provides a few more specialized quadrature rules. A list can be found on the [Quadrature Rules Namespace](@ref lf::quad) page.

# Using a Quadrature Rule

```cpp
auto ref_tria = lf::base::RefEl::kTria();

lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_tria, 3);

// Get degree of quadrature rule
int degree = qr.Degree();

// Get the order of a quadrature rule (degree + 1)
int order = qr.Order();

// Get points as columns of a matrix
Eigen::MatrixXd points = qr.Points();

// Get weights as a vector
Eigen::VectorXd weights = qr.Weights();
```

