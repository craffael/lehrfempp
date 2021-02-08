/**
 * @file
 * @brief Define the BrepSurface class
 * @author Raffael Casagrande
 * @date   2020-11-20 04:36:05
 * @copyright MIT License
 */

#ifndef __086252f2a1a544c48fc037db3f4fbb00
#define __086252f2a1a544c48fc037db3f4fbb00

#include "brep_geometry.h"

namespace lf::brep::interface {
class BrepSurface : public BrepGeometry {
 protected:
  BrepSurface() = default;
  BrepSurface(const BrepSurface&) = default;
  BrepSurface(BrepSurface&&) = default;
  BrepSurface& operator=(const BrepSurface&) = default;
  BrepSurface& operator=(BrepSurface&&) = default;

 public:
  ~BrepSurface() = default;

  [[nodiscard]] virtual Eigen::Vector3d GlobalSingle(
      const Eigen::Vector2d& local) const = 0;
  [[nodiscard]] virtual Eigen::Matrix<double, 3, 2> JacobianSingle(
      const Eigen::Vector2d& local) const = 0;
  [[nodiscard]] virtual std::pair<double, Eigen::Vector2d> Project(
      const Eigen::Vector3d& global) const = 0;
  [[nodiscard]] virtual bool IsInBoundingBox(
      const Eigen::Vector3d& global) const = 0;
  [[nodiscard]] virtual bool IsInside(const Eigen::Vector2d& local) const = 0;
};
}  // namespace lf::brep::interface

#endif  // __086252f2a1a544c48fc037db3f4fbb00
