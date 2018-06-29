
#include <gtest/gtest.h>
#include <lf/base/static_vars.h>

using namespace lf::base;

namespace lf::base::test {

class StaticVarsTest {
 public:
  StaticVarsTest() = default;
  ~StaticVarsTest(void) = default;

  static int ctrl_var_;
  static int other_var_;
};

CONTROLDECLARE(testvar, "testvar");
CLASSCONTROLDECLARE(StaticVarsTest, ctrl_var_, "ctrl_var");
CONTROLDECLARECOMMENT(StaticVarsTest, other_var_, "other_var", "A test case");

TEST(StaticVar, BasicTest) {
  ListCtrlVars(std::cout);
  SetCtrlVar("testvar", 42);
  StaticVarsTest::ctrl_var_ = 13;
  EXPECT_EQ(testvar, 42);
  std::cout << "After setting variables:" << std::endl;
  ListCtrlVars(std::cout);
}

TEST(StaticVar, FileTest) {
  SetCtrlVar("read_ctrl_vars_file", 1);
  ListCtrlVars(std::cout);
  ReadCtrlVarsFile("setup.vars");
  ListCtrlVars(std::cout);
}
}  // namespace lf::base::test
