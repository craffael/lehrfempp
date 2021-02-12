/**
 * @file
 * @brief Utility functions to unit test BrepGeometries.
 * @author Raffael Casagrande
 * @date   2021-01-29 10:13:02
 * @copyright MIT License
 */

#ifndef __bb9913668eeb4af49c42e857ac39e08e
#define __bb9913668eeb4af49c42e857ac39e08e

#include <gtest/gtest.h>
#include <lf/brep/interface/interface.h>

#include "lf/mesh/utils/all_codim_mesh_data_set.h"

namespace lf::brep::test_utils {

template <class F>
Eigen::Matrix3Xd approxJacobian(const F& f, Eigen::VectorXd local) {
  Eigen::Matrix3Xd result(3, local.rows());
  static const double eps = 1e-5;  // See Timothy Sauter, Numerical Analysis
  for (int i = 0; i < local.rows(); ++i) {
    Eigen::VectorXd dh =
        eps * Eigen::MatrixXd::Identity(local.rows(), local.rows()).col(i);
    result.col(i) = (f(local + dh) - f(local - dh)) / (2 * eps);
  }
  return result;
}

void CheckBrepGeometry(::lf::brep::interface::BrepGeometry const& geom,
                       Eigen::MatrixXd local);

}  // namespace lf::brep::test_utils

#endif  // __bb9913668eeb4af49c42e857ac39e08e
