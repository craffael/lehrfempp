
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
  static int arg_var_;
  static int output_ctrl_; // Added
};

CONTROLDECLARE(testvar, "testvar");
CLASSCONTROLDECLARE(StaticVarsTest, ctrl_var_, "ctrl_var");
CONTROLDECLARECOMMENT(StaticVarsTest, other_var_, "other_var", "A test case");
CONTROLDECLARECOMMENT(StaticVarsTest, arg_var_, "arg_var",
                      "Set from command line");

CONTROLDECLARECOMMENT(StaticVarsTest, output_ctrl_, "output_ctrl_", "Test output ctrl");

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
  SetCtrlVar("read_ctrl_vars_args", 10);
  std::cout << "########## Listing of control variables:" << std::endl;
  ListCtrlVars(std::cout);
  ReadCtrlVarsFile("setup.vars");
  std::cout << "########## After reading file:"
            << " Listing of control variables:" << std::endl;
  ListCtrlVars(std::cout);
  int argc = 2;
  const char* argv[] = {"arg_var=51", "another_option"};
  ReadCtrVarsCmdArgs(argc, argv);
  std::cout << "########## After reading command line args:"
            << " Listing of control variables:" << std::endl;
  ListCtrlVars(std::cout);
}

// Added to figure out how tests work
TEST(StaticVar, TestTest){
    ListCtrlVars(std::cout);
    SetCtrlVar("output_ctrl_", 5);
    StaticVarsTest::output_ctrl_ = 6;
    EXPECT_EQ(StaticVarsTest::output_ctrl_, 8)
            << "Output you get only if the test fails";
    std::cout << "After setting variables:" << std::endl;
    ListCtrlVars(std::cout);
}

}  // namespace lf::base::test
