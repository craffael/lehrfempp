

#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <utility>
#include "lf/base/random_access_iterator.h"

using namespace lf::base;

TEST(RandomAccessIterator, nonConst) {
  std::vector<int> numbers{0,1,2,3};
  RandomAccessIterator<int> ri = numbers.begin();

  EXPECT_EQ(*ri, 0);
  EXPECT_EQ(*++ri, 1);
  EXPECT_EQ(*ri, 1);
  EXPECT_EQ(*(ri++), 1);
  EXPECT_EQ(*ri, 2);
  EXPECT_EQ(ri, ri);

  RandomAccessIterator<int> ri2 = ri;
  EXPECT_EQ(ri, ri2);
  ++ri;
  EXPECT_NE(ri, ri2);
  ri = ri2;
  EXPECT_EQ(ri, ri2);
  EXPECT_EQ(*ri2, 2);
  EXPECT_EQ(*ri, 2);
  EXPECT_NE(ri==ri2, ri != ri2);

  auto ri3 = std::move(ri2);
  EXPECT_EQ(ri3, ri);

  *ri = 3;
  EXPECT_EQ( *ri3, 3);

  // the value of a const iterator can be modified.
  *static_cast<const RandomAccessIterator<int>>(ri) = 4;

  EXPECT_EQ(*ri, 4);
}

TEST(RandomAccessIterator, InteractWithForwardIterators) {
  std::vector<int> numbers{0,1,2,3};

  RandomAccessIterator<int> ri = numbers.begin();
  ForwardIterator<int> fi = ri;
  EXPECT_EQ(ri, fi);
  EXPECT_EQ(fi, ri);

  EXPECT_NE(ri += 2, fi);
  EXPECT_NE(ri, fi);
  EXPECT_EQ(*ri, 2);

  fi = std::move(ri);
  EXPECT_EQ(*fi, 2);
  ri = numbers.begin();
  EXPECT_EQ(*ri, 0);
  ri = ri + 2;
  EXPECT_EQ(*ri, 2);
  EXPECT_EQ(fi, ri);
}

TEST(RandomAccessIterator, RandomAccessFunctionality) {
  std::vector<int> numbers{0,1,2,3};
  RandomAccessIterator<int> ri0 = numbers.begin();
  auto ri1 = ri0 + 1;
  EXPECT_EQ(*ri1, 1);
  EXPECT_EQ(ri1-1, ri0);
  EXPECT_EQ(ri1-ri0, 1);
  EXPECT_EQ(ri0+1,ri1);
  EXPECT_TRUE(ri0 < ri1);
  EXPECT_FALSE(ri0 > ri1);
  EXPECT_TRUE(ri0 <= ri0);
  EXPECT_TRUE(ri1 >= ri1);

  ri0 += 1;
  EXPECT_EQ(ri0, ri1);
  ri0 -= 1;
  EXPECT_EQ(*ri0, 0);
}