/**
 * @file
 * @brief Implementation from occt_utils.h
 * @author Raffael Casagrande
 * @date   2020-11-11 02:31:37
 * @copyright MIT License
 */

#include "occt_details.h"

namespace lf::brep::occt::detail {
Eigen::Map<const Eigen::Vector3d> ToVector(const gp_Pnt& p) {
  return Eigen::Map<const Eigen::Vector3d>(
      reinterpret_cast<const double*>(&p.XYZ()));
}

Eigen::Map<const Eigen::Vector3d> ToVector(const gp_Vec& v) {
  return Eigen::Map<const Eigen::Vector3d>(
      reinterpret_cast<const double*>(&v.XYZ()));
}
}  // namespace lf::brep::occt::detail
