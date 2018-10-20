#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>

namespace lf::geometry::test {

TEST(QuadTest, checkJacobians) {
  // Define geometry of quadrilateral through corner positions
  Eigen::Matrix<double, Eigen::Dynamic, 4> quad_corners(2, 4);
  quad_corners << 2, 3, 3, 2, 2, 1, 2, 3;
  // Reference coordinates
  Eigen::MatrixXd refcoords(2, 2);
  refcoords << 0.2, 0.7, 0.1, 0.3;

  lf::geometry::QuadO1 generic_quad(quad_corners);
  lf::geometry::Parallelogram parallelogram(quad_corners);

  // Compute Jacobians
  const auto Jac_quad{generic_quad.Jacobian(refcoords)};
  const auto Jac_parg{parallelogram.Jacobian(refcoords)};
  EXPECT_NEAR((Jac_quad - Jac_parg).norm(), 0.0, 1e-12)
      << "Different Jacobians " << Jac_quad << " <-> " << Jac_parg;

  // Compute Inverse Gramians
  const auto Jacit_quad(generic_quad.JacobianInverseGramian(refcoords));
  const auto Jacit_parg(parallelogram.JacobianInverseGramian(refcoords));
  EXPECT_NEAR((Jacit_quad - Jacit_parg).norm(), 0.0, 1e-12)
      << "Different inverse Jacobians " << Jacit_quad << " <-> " << Jacit_parg;

  // Compare metric factors
  const auto g_quad{generic_quad.IntegrationElement(refcoords)};
  const auto g_parg{parallelogram.IntegrationElement(refcoords)};
  EXPECT_NEAR((g_quad - g_parg).norm(), 0.0, 1e-12)
      << "Different metric factors " << Jacit_quad << " <-> " << Jacit_parg;
}

}  // namespace lf::geometry::test
