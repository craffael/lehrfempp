#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <algorithm>

using namespace lf::base;

template<class T>
struct S {
  S(const std::list<T>& a) {}
};

TEST(Range, useInForLoop)
{
  std::list<int> list = {0,1,2,3};
  auto range = ForwardRange<int>(list.begin(),list.end());

  int count = 0;
  for(auto i : range) {
    EXPECT_EQ(i, count++);
  }


  // ForwardRange can be traversed twice
  count = 0;
  for(auto& i : range) {
    EXPECT_EQ(i, count++);
  }

  // try a few std:: algorithms
  EXPECT_EQ(*std::find(range.begin(), range.end(), 2), 2);
  EXPECT_EQ(std::count(range.begin(), range.end(), 2), 1);
  EXPECT_EQ(std::distance(range.begin(), range.end()), 4);
  bool result = std::all_of(range.begin(), range.end(), [](int a) {return true;});
  EXPECT_EQ(result, true);
}