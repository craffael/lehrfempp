#ifndef __b193e889f4e245959d124b0ded8a3b68
#define __b193e889f4e245959d124b0ded8a3b68

#include "geometry_interface.h"

namespace lf::geometry {

/** @brief Asserting a non-degenerate bilinear quadrilateral
 *
 * @param coords w x 4 Eigen matrix whose columns contain the vertex coordinates
 *        of the quadrilateral, w = world dimension
 * @param tol relative tolerance for numerical tests of equality with zero
 *
 * Terminates execution in case degenerate shape is detected. 
 */
bool assertNonDegenerateQuad(
    const Eigen::Matrix<double, Eigen::Dynamic, 4> &coords, double tol = 1.0E-8);

/**
 * @brief Bilinear quadrilateral element shape
 *
 * This element shape arises from subjecting the unit square to
 * a componentwise bilinear mapping. This results in a _non-flat_ quadrilateral
 * with _straight_ edges, whose shape is determined by the location of the
 * four vertices.
 */
class QuadO1 : public Geometry {
 public:
  /**
   * @brief Constructor building quadrilateral from vertex coordinates
   *
   * @param coords w x 4 matrix, w = world dimension, whose columns
   *        contain the world coordinates of the vertices.
   *
   * The constructor also checks whether the quadrilateral is non-degenerate
   */
  explicit QuadO1(Eigen::Matrix<double, Eigen::Dynamic, 4> coords);

  dim_t DimLocal() const override { return 2; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }

  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const ::Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;

  /** @copydoc Geometry::SubGeometry() */
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /** @copydoc Geometry::isAffine() */
  bool isAffine() const override { return false; }

  /**
   * @copydoc lf::geometry::Geometry::ChildGeometry()
   *
   * For a detailed description of the indexing of the vertices of child
   * entities see `Refinement.xoj`.
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat, lf::base::dim_t codim) const override;

 private:
  /** @brief Coordinates of the a four vertices, stored in matrix columns */
  Eigen::Matrix<double, Eigen::Dynamic, 4> coords_;
};

/**
 * @brief Affine quadrilateral = parallelogram
 *
 * Obtained as affine image of the unit square.
 */
class Parallelogram : public Geometry {
 public:
  /**
   * @brief Constructor building a parallelogram from _four_ vertex coordinates
   *
   * @param coords w x 4 matrix, w = world dimension, whose columns
   *        contain the world coordinates of the vertices.
   *
   * @note **third** vertex position is ignored, because it is already fixed
   *       by the three others in the case of a parallelogram
   */
  explicit Parallelogram(Eigen::Matrix<double, Eigen::Dynamic, 4> coords);
  /**
   * @brief Constructor building a parallelogram from _four_ vertex coordinates
   *
   * @param p0 coordinate vector for vertex 0
   * @param p1 coordinate vector for vertex 1
   * @param p2 coordinate vector for vertex 2
   *
   * Constructs a parallelogram from the coordinate vectors of three spanning
   * vertices.
   */
  explicit Parallelogram(const Eigen::VectorXd& p0, const Eigen::VectorXd& p1,
                         const Eigen::VectorXd& p2);

  dim_t DimLocal() const override { return 2; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }

  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const ::Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;

  /** @copydoc Geometry::SubGeometry() */
  std::unique_ptr<Geometry> SubGeometry(dim_t codim, dim_t i) const override;

  /** @copydoc Geometry::isAffine()
   *
   * @note In contrast to a general quadrilateral a parallelogram is
   *       is the **affine image** of a square.
   */
  bool isAffine() const override { return true; }

  /**
   * @copydoc lf::geometry::Geometry::ChildGeometry()
   *
   * For a detailed description of the indexing of the vertices of child
   * entities see `Refinement.xoj`.
   */
  std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pat, lf::base::dim_t codim) const override;

 private:
  /** @brief performs initialization of data members */
  void init(void);

  /** @brief Coordinates of the a four vertices, stored in matrix columns
   *
   * In fact, three vertex coordinates of the _spanning vertices_ would be
   * sufficient to define the shape of a parallelogram. For convenience all
   * vertices are stored as in QuadO1.
   */
  Eigen::Matrix<double, Eigen::Dynamic, 4> coords_;
  /** @brief Matrix of affine mapping generating the parallelogram,
   *        agrees with constant Jacobian */
  Eigen::Matrix<double, Eigen::Dynamic, 2> jacobian_;
  /** @brief Constant matrix (\f$ J (J^T J)^{-1}\f$) */
  Eigen::Matrix<double, Eigen::Dynamic, 2> jacobian_inverse_gramian_;
  /** @brief constant Gramian determinant */
  double integrationElement_;
};

}  // namespace lf::geometry

#endif  // __b193e889f4e245959d124b0ded8a3b68
