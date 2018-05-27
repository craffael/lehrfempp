
#ifndef __7ed6b0d4d9244155819c464fc4eb9bbb
#define __7ed6b0d4d9244155819c464fc4eb9bbb

#include <lf/base/ref_el.h>
#include <Eigen/Eigen>

namespace lf::geometry {

class Geometry {
 public:


  /**
   * @brief Dimension of the domain of this mapping.
   */
  virtual char DimLocal() const = 0;


  /**
   * @brief Dimension of the image of this mapping.
   */
  virtual char DimGlobal() const = 0;


  /**
   * @brief The Reference element that defines the domain of this mapping.
   */
  virtual base::RefEl RefEl() const = 0;


  /**
   * @brief Map a number of points in local coordinates into the global 
   *        coordinate system.
   * @param local A Matrix of size `DimFrom() x numPoints` that contains
   *              the points at which the mapping function should be evaluated
   *              as column vectors.
   * @return A Matrix of size `DimTo() x numPoints` that contains the mapped
   *         points as column vectors.
   */
  virtual Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const = 0;


  /**
   * @brief Evaluate the jacobian of the mapping simultaneously at `numPoints`
   *        points.
   * @param local A Matrix of size `DimFrom() x numPoints` that contains the 
   *              evaluation points as column vectors
   * @return A Matrix of size `DimTo() x (DimFrom() * numPoints)` that contains
   *         the jacobians at the evaluation points.
   */
  virtual Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const = 0;


  /**
   * @brief The integration element (factor appearing in integral transformation
   *        formula, see below) at number of evaluation points (specified in
   *        local coordinates).
   * @param local A Matrix of size `DimFrom() x numPoints` that contains the 
   *              evaluation points (in local coordinates) as column vectors.
   * @return A Vector of size `numPoints x 1` that contains the integration
   *         elements at every evaluation point.
   *         
   * For a transformation \f$ \Phi : K \mapsto R^{\text{dimTo}}\f$ with Jacobian 
   * \f$ D\Phi : K \mapsto R^{\text{dimTo} \times \text{dimFrom}} \f$ the integration
   * element \f$ g \f$ at point \f$ \xi \in K \f$ is defined by
   * \f[
   * g(\xi) := \sqrt{\mathrm{det}\left|D\Phi^T(\xi) D\Phi(\xi) \right|}
   * \f]
   */
  virtual Eigen::VectorXd IntegrationElement(const Eigen::MatrixXd& local) const = 0;


  /**
   * @brief Virtual destructor
   */
  virtual ~Geometry() = default;
};

}  // namespace lf::mesh

#endif  // __7ed6b0d4d9244155819c464fc4eb9bbb
