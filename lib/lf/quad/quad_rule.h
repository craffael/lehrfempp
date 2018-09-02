/**
 * @file
 * @brief Declaration of the QuadRule class
 * @author Raffael Casagrande
 * @date   2018-08-11 03:25:35
 * @copyright MIT License
 */

#ifndef __a7241ee797424d98ad339341b02bca70
#define __a7241ee797424d98ad339341b02bca70

#include <lf/base/base.h>
#include <utility>

namespace lf::quad {

using quadOrder_t = unsigned char;

/**
 * @brief Represents a Quadrature Rule over one of the Reference Elements
 *
 * A Quadrature rule essentially consists of quadrature nodes \f$ \{\vec{\xi}_0,
 * \ldots, \vec{\xi}_{n-1}\} \f$ and quadrature weights \f$ \{\omega_1, \ldots,
 * \omega_{n-1}\}\f$.
 *
 * The integral of a function \f$ f \f$ over the reference element \f$K\f$
 * is then approximated by
 * \f[
 * \int_K f(\vec{x}) \, d\vec{x} \approx \sum_{i=0}^{n-1} f(\vec{\xi_i})
 * \omega_i
 * \f]
 *
 */
class QuadRule {
 public:
  /**
   * @brief Construct a new quadrature rule by specifying reference element,
   * points, weights and order explicitly.
   * @param ref_el The reference element for which the quadrature rule is.
   * @param points The points of the quadrature rule, a matrix of size
   * `ref_el.Dimension() x num_points` that contains the points as column
   * vectors
   * @param weights The weights of the quadrature rule, a vector of length
   * `num_points`
   * @param order The order of the quadrature rule, see Order() for more info.
   */
  explicit QuadRule(base::RefEl ref_el, Eigen::MatrixXd points,
                    Eigen::VectorXd weights, quadOrder_t order)
      : ref_el_(std::move(ref_el)),
        order_(order),
        points_(std::move(points)),
        weights_(std::move(weights)) {
    LF_ASSERT_MSG(points_.rows() == ref_el_.Dimension(), "Dimension mismatch");
    LF_ASSERT_MSG(weights_.size() == points_.cols(), "Dimension mismatch");
  }

  /**
   * @brief The reference element \f$ K \f$ over which this QuadRule
   * integrates.
   */
  base::RefEl RefEl() const { return ref_el_; }

  /**
   * @brief Return the order of this Quadrature Rule.
   *
   * The order of a quadrature rule is defined as the largest integer \f$k\f$
   * such that \f[
   * \forall p \in \mathbb{P}_k, \quad \int_K p(\vec{x}) \,
   * d\vec{x} = \sum_{i=0}^{n-1} \omega_i f(\vec{\xi_i})
   * \f]
   * here the polynomial space \f$\mathbb{P}_k\f$ is defined differently for
   * every reference element:
   *
   * - For base::RefEl::kSegment(), \f$ \mathbb{P}_k := \mathrm{span} \{ x^a
   * \, | \, 0 \leq a \leq k \} \f$.
   * - For base::RefEl::kTria(), \f$ \mathbb{P}_k := \mathrm{span} \{ x^a y^b \,
   * | \, 0 \leq a + b \leq k \} \f$
   * - For base::RefEl::kQuad(), \f$ \mathbb{P}_k := \mathrm{span} \{ x^a y^b \,
   * | \, 0 \leq \mathrm{max}(a,b) \leq k \} \f$
   *
   */
  quadOrder_t Order() const { return order_; }

  /**
   * @brief All quadrature points \f$ \begin{pmatrix} \vec{\xi}_0, \ldots,
   * \vec{\xi}_{n-1}  \end{pmatrix} \f$ as column vectors.
   * @return A matrix of size `RefEl().Dimension() x ` \f$ n \f$ that contains
   *         the point coordinates as column vectors.
   */
  const Eigen::MatrixXd& Points() const { return points_; }

  /**
   * @brief All quadrature weights \f$ \begin{pmatrix} \omega_0, \ldots,
   * \omega_{n-1} \end{pmatrix} \f$ as a vector.
   * @return A vector of length \f$ n \f$
   */
  const Eigen::VectorXd& Weights() const { return weights_; }

  /**
   * @brief Return the total number \f$ n \f$ of quadrature points (num of
   * columns of points/weights)
   */
  base::size_type NumPoints() const { return points_.cols(); }

 private:
  base::RefEl ref_el_;
  quadOrder_t order_;
  Eigen::MatrixXd points_;
  Eigen::VectorXd weights_;
};

}  // namespace lf::quad

#endif  // __a7241ee797424d98ad339341b02bca70
