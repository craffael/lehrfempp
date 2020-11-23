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
#include "brep_surface.h"
#include "lf/base/base.h"
#include "lf/mesh/utils/all_codim_mesh_data_set.h"

namespace lf::brep::interface {

class BrepModel {
 protected:
  BrepModel() = default;

 public:
  BrepModel(const BrepModel&) = delete;
  BrepModel(BrepModel&&) = delete;
  BrepModel& operator=(const BrepModel&) = delete;
  BrepModel& operator=(BrepModel&&) = delete;
  virtual ~BrepModel() = default;

  [[nodiscard]] virtual std::vector<std::pair<BrepCurve const*, double>>
  FindCurves(const Eigen::Vector3d& global) const = 0;

  [[nodiscard]] virtual std::vector<
      std::pair<BrepCurve const*, Eigen::RowVectorXd>>
  FindCurvesMulti(const Eigen::Matrix3Xd& global) const = 0;

  [[nodiscard]] virtual std::vector<
      std::pair<BrepSurface const*, Eigen::Vector2d>>
  FindSurfaces(const Eigen::Vector3d& global) const = 0;

  [[nodiscard]] virtual std::vector<
      std::pair<BrepSurface const*, Eigen::Matrix2Xd>>
  FindSurfacesMulti(const Eigen::Matrix3Xd& global) const = 0;

  [[nodiscard]] virtual base::size_type NumCurves() const = 0;
  [[nodiscard]] virtual base::size_type NumSurfaces() const = 0;
};

}  // namespace lf::brep::interface

#endif  // __c7d5be9b848449af86de926a04d5d43c
