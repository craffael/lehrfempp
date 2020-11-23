/**
 * @file
 * @brief Defines the class BrepCurve
 * @author Raffael Casagrande
 * @date   2020-11-20 05:16:31
 * @copyright MIT License
 */

#ifndef __3e510e104ccc4effa7d2989488449f5a
#define __3e510e104ccc4effa7d2989488449f5a
#include "brep_geometry.h"

namespace lf::brep::interface {

class BrepCurve : public BrepGeometry {
 protected:
  BrepCurve() = default;
  BrepCurve(const BrepCurve&) = default;
  BrepCurve(BrepCurve&&) = default;
  BrepCurve& operator=(const BrepCurve&) = default;
  BrepCurve& operator=(BrepCurve&&) = default;

 public:
  ~BrepCurve() = default;

  [[nodiscard]] virtual Eigen::Vector3d GlobalSingle(double local) const = 0;
  [[nodiscard]] virtual Eigen::Vector3d JacobianSingle(double local) const = 0;

  
  /**
   * @brief Project a point onto this curve and retrieve the local coordinates.
   * @param global The (global) coordinates of the point that should be projeted.
   * @return The distance of the projected point to the original point (first)
   *         followed by the local coordinate of the point w.r.t. the curve.
   *
   * It's important to realize, that Project may return a point that lies on the underlying
   * curve but that is "outside" the parameter bounds, i.e. `IsInside() == false`.
   *
   * If `IsInside(local) == true`, we have the following identity (up to rounding errors):
   * `Project(Global(local)) == local`.
   *
   * @note This does not necessarily hold if `IsInside(local)==false`!
   *
   *
   */
  [[nodiscard]] virtual std::pair<double, double> Project(
      const Eigen::Vector3d& global) const = 0;
  [[nodiscard]] virtual bool IsInBoundingBoxSingle(
      const Eigen::Vector3d& global) const = 0;
  [[nodiscard]] virtual bool IsInside(double local) const = 0;
};

}  // namespace lf::brep::interface

#endif  // __3e510e104ccc4effa7d2989488449f5a
