/**
 * @file
 * @brief Some utility functions to deal with OpenCascade, should not be used externally.
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
  LF_ASSERT_MSG(m.innerSize() == 3, "m must be 3d vector.");
  return gp_Pnt(m[0], m[1], m[2]);
}

template <class DERIVED>
gp_Pnt2d ToPoint2d(const Eigen::DenseBase<DERIVED>& m) {
  LF_ASSERT_MSG(m.outerSize() == 1, "m must be a vector");
  LF_ASSERT_MSG(m.innerSize() == 2, "m must be a 2d vector");
  return gp_Pnt2d(m[0], m[1]);
}

Eigen::Map<const Eigen::Vector3d> ToVector(const gp_Pnt& p);

Eigen::Map<const Eigen::Vector3d> ToVector(const gp_Vec& v);

}  // namespace lf::brep::occt::detail

#endif  // __dd390c8ad3094cf296921339f483f15b
