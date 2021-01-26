/**
 * @file
 * @brief Test implementation of QuadRuleCache
 * @author Raffael Casagrande
 * @date   2021-01-14 02:46:42
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/quad/quad.h>

namespace lf::quad::test {

TEST(lf_quad_QuadRuleCache, Cache) {
  QuadRuleCache qrc;
  for (auto ref_el :
       {base::RefEl::kSegment(), base::RefEl::kTria(), base::RefEl::kQuad()}) {
    for (int p = 0; p < 10; ++p) {
      auto& qr1 = qrc.Get(ref_el, p);
      auto qr1_ref = make_QuadRule(ref_el, p);
      EXPECT_EQ(qr1.NumPoints(), qr1_ref.NumPoints());
      EXPECT_EQ(qr1.RefEl(), ref_el);
      EXPECT_EQ(qr1.Degree(), qr1_ref.Degree());
      EXPECT_TRUE(qr1.Points().isApprox(qr1_ref.Points()));
      EXPECT_TRUE(qr1.Weights().isApprox(qr1_ref.Weights()));

      auto& qr11 = qrc.Get(ref_el, p);
      EXPECT_EQ(&qr11, &qr1);
    }
  }
}

TEST(lf_quad_QuadRuleCache, ReferencesStayValid) {
  QuadRuleCache qrc;
  const auto& qr1 = qrc.Get(base::RefEl::kSegment(), 3);
  auto qr1_val = qr1;
  // ask for a few more quadrature rules:
  for (int p = 3; p < 10; ++p) {
    auto q = qrc.Get(base::RefEl::kSegment(), p);
  }
  // make sure the old reference is still valid:
  auto& qr1_new = qrc.Get(base::RefEl::kSegment(), 3);
  EXPECT_EQ(&qr1, &qr1_new);
}

}  // namespace lf::quad::test
