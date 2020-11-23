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
  [[nodiscard]] virtual std::pair<double, double> Project(
      const Eigen::Vector3d& global) const = 0;
  [[nodiscard]] virtual bool IsInBoundingBoxSingle(
      const Eigen::Vector3d& global) const = 0;
  [[nodiscard]] virtual bool IsInside(double local) const = 0;
};

}  // namespace lf::brep::interface

#endif  // __3e510e104ccc4effa7d2989488449f5a
