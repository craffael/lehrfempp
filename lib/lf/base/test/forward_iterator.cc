#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <utility>

using namespace lf::base;

TEST(ForwardIterator, nonConst) {
  std::vector<int> numbers{0, 1, 2, 3};
  ForwardIterator<int> fi = numbers.begin();
  EXPECT_EQ(*fi, 0);
  EXPECT_EQ(*++fi, 1);
  EXPECT_EQ(*fi, 1);
  EXPECT_EQ(*(fi++), 1);
  EXPECT_EQ(*fi, 2);
  EXPECT_EQ(fi, fi);

  ForwardIterator<int> fi2 = fi;
  EXPECT_EQ(fi, fi2);
  ++fi;
  EXPECT_NE(fi, fi2);
  fi = fi2;
  EXPECT_EQ(fi, fi2);
  EXPECT_EQ(*fi2, 2);
  EXPECT_EQ(*fi, 2);
  EXPECT_NE(fi==fi2, fi != fi2);

  auto fi3 = std::move(fi2);
  EXPECT_EQ(fi3, fi);

  *fi = 3;
  EXPECT_EQ( *fi3, 3);
  *std::as_const(fi) = 4; // the value of a const iterator can be modified.
  EXPECT_EQ(*fi, 4);

  ForwardIterator<const int> fiConst = fi;
  EXPECT_EQ(*fiConst, 4);
}

TEST(ForwardIterator, constEntities) {
  const std::vector<std::string> strings{"hello", "world"};

  ForwardIterator<const std::string> fi = strings.begin();
  EXPECT_EQ(*fi, "hello");
  EXPECT_EQ(*fi++, "hello");
  EXPECT_EQ(*fi, "world");
  ++fi;
  EXPECT_EQ(fi, strings.end());
}

TEST(ForwardIterator, DefaultConstructible) {
  std::vector<int> numbers{0,1,2};

  ForwardIterator<int> fi0 = numbers.begin();
  ForwardIterator<int> defaultFi;
  EXPECT_NE(fi0, defaultFi);
  fi0 = defaultFi;
  EXPECT_EQ(defaultFi, fi0);

  defaultFi = numbers.end();
  EXPECT_NE(defaultFi, fi0);
}
