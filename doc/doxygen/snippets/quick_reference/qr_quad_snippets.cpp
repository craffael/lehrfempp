// This file contains Doxygen snippets for the quick reference quad document
// It defines a function to hold all code snippets and includes necessary
// imports.

#include <lf/base/base.h>
#include <lf/quad/quad.h>

#include <iostream>

void qr_geometry_snippets() {
  {
    //! [quad Working with Quadrature Rules]
    lf::base::RefEl ref_tria = lf::base::RefEl::kTria();

    lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_tria, 3);
    //! [quad Working with Quadrature Rules]
  }

  {
    //! [quad Working with Quadrature Rules 0]
    lf::base::RefEl ref_tria = lf::base::RefEl::kTria();

    lf::quad::QuadRule qr = lf::quad::make_QuadRuleNodal(ref_tria);
    //! [quad Working with Quadrature Rules 0]
  }

  {
    //! [quad Working with Quadrature Rules 1]
    lf::base::RefEl ref_tria = lf::base::RefEl::kTria();

    Eigen::MatrixXd points(2, 3);
    points << 0.166667, 0.666667, 0.166667, 0.166667, 0.166667, 0.666667;

    Eigen::Vector3d weights;
    weights << 0.166667, 0.166667, 0.166667;

    lf::quad::QuadRule qr(ref_tria, points, weights, 2);
    //! [quad Working with Quadrature Rules 1]
  }

  {
    //! [quad Using a Quadrature Rule]
    lf::base::RefEl ref_tria = lf::base::RefEl::kTria();

    lf::quad::QuadRule qr = lf::quad::make_QuadRule(ref_tria, 3);

    // Get degree of quadrature rule
    int degree = qr.Degree();

    // Get the order of a quadrature rule (degree + 1)
    int order = qr.Order();

    // Get points as columns of a matrix
    Eigen::MatrixXd points = qr.Points();

    // Get weights as a vector
    Eigen::VectorXd weights = qr.Weights();
    //! [quad Using a Quadrature Rule]
  }
}
