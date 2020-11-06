/**
 * @file
 * @brief Describes the BrepModel abstract interface class
 * @author Raffael Casagrande
 * @date   2020-11-06 04:43:15
 * @copyright MIT License
 */

#ifndef __c7d5be9b848449af86de926a04d5d43c
#define __c7d5be9b848449af86de926a04d5d43c

#include "brep_geometry.h"
#include "lf/base/base.h"

namespace lf::brep {

class BrepModel {
 protected:
  BrepModel() = default;

 public:
  BrepModel(const BrepModel&) = delete;
  BrepModel(BrepModel&&) = delete;
  BrepModel& operator=(const BrepModel&) = delete;
  BrepModel& operator=(BrepModel&&) = delete;
  virtual ~BrepModel() = default;

  [[nodiscard]] virtual std::pair<std::shared_ptr<BrepGeometry>,
                                  Eigen::MatrixXd>
  FindGeometry(const Eigen::MatrixXd& global) const = 0;
};

}  // namespace lf::brep

#endif  // __c7d5be9b848449af86de926a04d5d43c
