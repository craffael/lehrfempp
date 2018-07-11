

#include <gtest/gtest.h>
#include <lf/base/base.h>

using namespace lf::base;

namespace lf::base::test {
TEST(RefEl, dimensionCorrect) {
  EXPECT_EQ(RefEl::kPoint().Dimension(), 0);
  EXPECT_EQ(RefEl::kSegment().Dimension(), 1);
  EXPECT_EQ(RefEl::kTria().Dimension(), 2);
  EXPECT_EQ(RefEl::kQuad().Dimension(), 2);
}

TEST(RefEl, numberOfNodesCorrect){
    EXPECT_EQ(RefEl::kPoint().NumNodes(), 1);
    EXPECT_EQ(RefEl::kSegment().NumNodes(), 2);
    EXPECT_EQ(RefEl::kTria().NumNodes(), 3);
    EXPECT_EQ(RefEl::kQuad().NumNodes(), 4);

}

}  // namespace lf::base::test
