/**
 * @file
 * @brief Data structures representing HP finite elements
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include "hierarchic_fe.h"

namespace lf::fe {

Eigen::VectorXd chebyshevNodes(unsigned n) {
  // Compute the chebyshev nodes in the interval [0, 1]
  const auto cosine = [](double x) -> double { return std::cos(x); };
  return (Eigen::VectorXd::Ones(n) +
          Eigen::VectorXd::LinSpaced(n, M_PI - M_PI / (2 * n), M_PI / (2 * n))
              .unaryExpr(cosine)) /
         2;
}

}  // end namespace lf::fe
