/**
 * @file
 * @brief Implementations from make_quad_rule.h
 * @author Raffael Casagrande
 * @date   2018-08-19 06:35:16
 * @copyright MIT License
 */

#include "make_quad_rule.h"
#include <Eigen/KroneckerProduct>
#include "gauss_legendre.h"

namespace lf::quad {
QuadRule make_QuadRule(base::RefEl ref_el, unsigned char order) {
  if (ref_el == base::RefEl::kSegment()) {
    unsigned int n = order / 2 + 1;
    auto [points, weights] = GaussLegendre(n);
    return QuadRule{base::RefEl::kSegment(), std::move(points),
                    std::move(weights), 2 * n - 1};
  }
  if (ref_el == base::RefEl::kQuad()) {
    unsigned int n = order / 2 + 1;
    auto [points1d, weights1d] = GaussLegendre(n);
    Eigen::MatrixXd points2d(n * n);
    points2d.row(0) =
        Eigen::kroneckerProduct(points1d, Eigen::VectorXd::Ones(n));
    points2d.row(1) = points1d.replicate(n);
    return QuadRule{base::RefEl::kQuad(), std::move(points2d),
                    Eigen::kroneckerProduct(weights1d, weights1d), 2 * n - 1};
  }
  LF_VERIFY_MSG(
      false, "No Quadrature rules implemented for this reference element yet.");
}
}  // namespace lf::quad
