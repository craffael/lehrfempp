## Solution of General Second-Order Elliptic BVP by means of linear Finite Elements

File: `ellbvp_linfe_demo.cc`

This demonstration code
- relies on a triangular computational domain,
- defines all the coefficient functions and source functions for a rather general
  second-order elliptic boundary value problem,
- builds a simple 2D hybrid mesh comprising triangles and quadrilaterals,
- constructs a hierarchy of nested meshes by uniform refinement,
- assembles the Galerkin matrix and right hand side vector on each refinement level
  using lowest-order (linear) Lagrangian finite elements,
- solves the finite element linear system,
- approximately computes the L2 norm and H1 seminorm of the finite element discretization
  errors

File: `homDir_linfe_demo.cc`

This small demo shows how to solve a second-order linear scalar boundary value problem
with zero Dirichlet boundary conditions by means a Galerkin finite-element scheme based on
piecewise linear Lagrangian finite elements. 
