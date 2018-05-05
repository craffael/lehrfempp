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
  auto range = Range<int>(list.begin(),list.end());

  int count = 0;
  for(auto i : range) {
    EXPECT_EQ(i, count++);
  }


  // Range can be traversed twice
  count = 0;
  for(auto& i : range) {
    EXPECT_EQ(i, count++);
  }

  // std:: algorithms work
  //EXPECT_EQ(std::find(range.begin(), range.end(), 2), 2);

}