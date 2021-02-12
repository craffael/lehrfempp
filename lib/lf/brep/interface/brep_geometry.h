/**
 * @file
 * @brief Contains the BrepGeometry interface class
 * @author Raffael Casagrande
 * @date   2020-11-06 04:39:39
 * @copyright MIT License
 */

#ifndef __b67c777b4ad848da8a7a5125cbff9b73
#define __b67c777b4ad848da8a7a5125cbff9b73

#include "lf/base/base.h"

namespace lf::brep::interface {

class BrepGeometry {
 protected:
  BrepGeometry() = default;

  BrepGeometry(const BrepGeometry&) = default;
  BrepGeometry(BrepGeometry&&) = default;
  BrepGeometry& operator=(const BrepGeometry&) = default;
  BrepGeometry& operator=(BrepGeometry&&) = default;

 public:
  virtual ~BrepGeometry() = default;

  [[nodiscard]] virtual base::dim_t DimGlobal() const = 0;
  [[nodiscard]] virtual base::dim_t DimLocal() const = 0;

  [[nodiscard]] virtual Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const = 0;
  [[nodiscard]] virtual Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const = 0;
  [[nodiscard]] virtual std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const = 0;

  /**
   * @brief Project a point onto this curve and retrieve the local coordinates.
   * @param global The (global) coordinates of the point that should be
   * projeted.
   * @return The distance of the projected point to the original point (first)
   *         followed by the local coordinate of the point w.r.t. the curve.
   *
   * It's important to realize, that Project may return a point that lies on the
   * underlying curve but that is "outside" the parameter bounds, i.e.
   * `IsInside() == false`.
   *
   * If `IsInside(local) == true`, we have the following identity (up to
   * rounding errors): `Project(Global(local)) == local`.
   *
   * @note This does not necessarily hold if `IsInside(local)==false`!
   *
   *
   */
  [[nodiscard]] virtual std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd& global) const = 0;

  /**
   * @brief Check if the given local coordinate is "inside"
   *
   * Note that the parametrization of the curve may be defined for all Real
   * numbers. But most of the times, the actual CurveEntity is finite and thus
   * the parameters are bounded. You can check whether the parameter is inside
   * the bounds using this function.
   *
   * @note This method uses a tolerance value to make sure that it also return
   * `true` for the endpoints of the curve.
   *
   * @note For periodic curves, such as a circle, this method will assume
   * arbitrary parameter boundaries that span a period (such as [0, 2*pi])
   */
  [[nodiscard]] virtual bool IsInside(const Eigen::VectorXd& local) const = 0;

  /**
   * @brief Return periods (in local coordinates) of this geometry. If the
   * geometry is not periodic at all, returns a vector with all elements exactly
   * zero.
   * @return A vector of length `DimLocal()` that contains the periods for every
   * axis.
   *
   * A Geometry is periodic with period `T` along axis `e_i` if `Global(x)` is
   * the same (up to rounding errors) as `Global(x+T*e_i)`
   */
  [[nodiscard]] virtual Eigen::VectorXd Periods() const = 0;
};

}  // namespace lf::brep::interface

#endif  // __b67c777b4ad848da8a7a5125cbff9b73
