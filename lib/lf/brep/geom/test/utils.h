/**
 * @file
 * @brief test utilities for lf::brep::geom::test
 * @author Raffael Casagrande
 * @date   2021-02-01 04:08:02
 * @copyright MIT License
 */

#ifndef __0be58301f9e84b5c8e76f46dc22888e5
#define __0be58301f9e84b5c8e76f46dc22888e5

#include <lf/geometry/refinement_pattern.h>

namespace lf::brep::geom::test {

/**
 * @brief Divides the Segment into two child segments
 */
class SegmentRefPat : public geometry::RefinementPattern {
 public:
  SegmentRefPat() : geometry::RefinementPattern(base::RefEl::kSegment(), 12) {}
  [[nodiscard]] lf::base::size_type NumChildren(
      lf::base::dim_t codim) const override {
    LF_ASSERT_MSG(codim < 2, "codim out of range.");
    return codim == 0 ? 2 : 3;
  }

  [[nodiscard]] std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
  ChildPolygons(lf::base::dim_t codim) const override {
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> result;
    if (codim == 0) {
      result.resize(2, Eigen::MatrixXi(1, 2));
      result[0] << 0, 6;
      result[1] << 6, 12;
    } else {
      result.resize(3, Eigen::MatrixXi(1, 1));
      result[0] << 0;
      result[1] << 6;
      result[2] << 12;
    }
    return result;
  }
};

}  // namespace lf::brep::geom::test

#endif  // __0be58301f9e84b5c8e76f46dc22888e5
