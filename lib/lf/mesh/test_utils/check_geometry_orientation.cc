#include "check_geometry_orientation.h"
#include <gtest/gtest.h>

namespace lf::mesh::test_utils {
void checkGeometryOrientation(const Entity& e) {
  using dim_t = base::RefEl::dim_t;
  auto zero = Eigen::VectorXd::Zero(0);
  for (dim_t i = 0; i < e.RefEl().NumNodes(); ++i) {
    base::RefEl ref_el = e.RefEl();
    Eigen::VectorXd coorda = e.Geometry()->Global(ref_el.NodeCoords().col(i));
    Eigen::VectorXd coordb =
        e.SubEntities(ref_el.Dimension())[i]->Geometry()->Global(zero);
    for (dim_t j = 0; j < coorda.rows(); ++j) {
      EXPECT_EQ(coorda(j), coordb(j)) << "a(" << j << ") != b(" << j << ")";
    }
  }
}
}  // namespace lf::mesh::test_utils
