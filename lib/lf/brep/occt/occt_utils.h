/**
 * @file
 * @brief Some utility functions to deal with OpenCascade
 * @author Raffael Casagrande
 * @date   2020-11-11 02:28:28
 * @copyright MIT License
 */

#ifndef __dd390c8ad3094cf296921339f483f15b
#define __dd390c8ad3094cf296921339f483f15b
#include <Eigen/Core>
#include <gp_Pnt.hxx>

#include "lf/base/base.h"

namespace lf::brep::occt::detail {
template <class DERIVED>
gp_Pnt ToPoint(const Eigen::DenseBase<DERIVED>& m) {
  LF_ASSERT_MSG(m.outerSize() == 1, "m must be a vector.");
  LF_ASSERT_MSG(m.innerSize() == 2 || m.innerSize() == 3,
                "m must be 2d or 3d vector.");
  gp_Pnt p;
  if (m.rows() == 2) {
    p = gp_Pnt(m[0], m[1], 0.0);
  } else {
    p = gp_Pnt(m[0], m[1], m[2]);
  }
  return p;
}

Eigen::Map<const Eigen::Vector3d> ToVector(const gp_Pnt& p);

Eigen::Map<const Eigen::Vector3d> ToVector(const gp_Vec& v);

}  // namespace lf::brep::occt::detail

#endif  // __dd390c8ad3094cf296921339f483f15b
