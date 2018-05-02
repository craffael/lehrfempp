


#include <gtest/gtest.h>
#include <lf/base/base.h>

using namespace lf::base;

TEST(RefEl, dimensionCorrect)
{
  EXPECT_EQ(RefEl::kPoint.Dimension(), 0);
  EXPECT_EQ(RefEl::kSegment.Dimension(), 1);
  EXPECT_EQ(RefEl::kTria.Dimension(), 2);
  EXPECT_EQ(RefEl::kQuad.Dimension(), 3);
}