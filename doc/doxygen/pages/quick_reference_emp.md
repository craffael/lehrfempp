# Quick Reference - Enitiy Matrix Providers {#quick_reference_emp}

[TOC]

> [!caution]
> Discussed in [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf) @lref{rem:cppconcepts}. Please read before using quick reference.

## Overview

<!-- Quick explanation and reference to relevant lecture document chapter -->
LehrFEM++ allows users to implement custom

## Custom Entity Vector Providers

```cpp

```

## Built-in Entity Matrix Providers

<!-- TODO (barmstron): Write down quick specification of each EL/VE MP -->

### lf::fe::DiffusionElementMatrixProvider
### lf::fe::MassElementMatrixProvider
### lf::fe::MassEdgeMatrixProvider
### lf::uscalfe::LinearFELaplaceElementMatrix
### lf::uscalfe::ReactionDiffusionElementMatrixProvider
### lf::uscalfe::MassEdgeMatrixProvider

## Built-in Entity Vector Providers

- lf::fe::ScalarLoadElementVectorProvider
- lf::fe::ScalarLoadEdgeVectorProvider
- lf::uscalfe::LinearFELocalLoadVector
- lf::uscalfe::ScalarLoadElementVectorProvider
- lf::uscalfe::ScalarLoadEdgeVectorProvider
