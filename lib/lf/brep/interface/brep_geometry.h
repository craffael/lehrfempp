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

namespace lf::brep {

class BrepGeometry {
 protected:
  BrepGeometry() = default;

 public:
  BrepGeometry(const BrepGeometry&) = delete;
  BrepGeometry(BrepGeometry&&) = delete;
  BrepGeometry& operator=(const BrepGeometry&) = delete;
  BrepGeometry& operator=(BrepGeometry&&) = delete;
  virtual ~BrepGeometry() = default;

  [[nodiscard]] virtual base::dim_t DimGlobal() const = 0;
  [[nodiscard]] virtual base::dim_t DimLocal() const = 0;
  [[nodiscard]] virtual Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const = 0;
  [[nodiscard]] virtual Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const = 0;
  [[nodiscard]] virtual std::pair<Eigen::VectorXd, Eigen::MatrixXd> Project(
      const Eigen::MatrixXd& local) const = 0;
  [[nodiscard]] virtual std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const = 0;
};

}  // namespace lf::brep

#endif  // __b67c777b4ad848da8a7a5125cbff9b73
