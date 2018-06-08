#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <utility>

namespace lf::base::test {

TEST(ForwardIteratorTest, nonConst) {
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
  EXPECT_NE(fi == fi2, fi != fi2);

  auto fi3 = std::move(fi2);
  EXPECT_EQ(fi3, fi);

  *fi = 3;
  EXPECT_EQ(*fi3, 3);

  // the value of a const iterator can be modified.
  *static_cast<const ForwardIterator<int>>(fi) = 4;

  EXPECT_EQ(*fi, 4);
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
  std::vector<int> numbers{0, 1, 2};

  ForwardIterator<int> fi0 = numbers.begin();
  ForwardIterator<int> defaultFi;
  EXPECT_NE(fi0, defaultFi);
  fi0 = defaultFi;
  EXPECT_EQ(defaultFi, fi0);

  defaultFi = numbers.end();
  EXPECT_NE(defaultFi, fi0);
}
}  // namespace lf::base::test
