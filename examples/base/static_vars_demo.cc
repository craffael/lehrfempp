/** @file static_vars_demo.cc
 * @brief Short demo showing how to set and use static variables
 */

#include <iostream>
#include "lf/base/base.h"
#include "lf/geometry/geometry.h"
#include "lf/mesh/mesh.h"
#include "lf/refinement/refinement.h"

class StaticVarsDemoClass {
 public:
  static int ctrl_var_;
  static int other_var_;
  static int arg_var_;
  static int output_ctrl_;  // Added
};

// Definition of static control variables

CONTROLDECLARE(testvar, "testvar");

CLASSCONTROLDECLARE(StaticVarsDemoClass, ctrl_var_, "ctrl_var");

CONTROLDECLARECOMMENT(StaticVarsDemoClass, other_var_, "other_var",
                      "A test case");

CONTROLDECLARECOMMENT(StaticVarsDemoClass, arg_var_, "arg_var",
                      "Set from command line");

CONTROLDECLARECOMMENT(StaticVarsDemoClass, output_ctrl_, "output_ctrl_",
                      "Test output ctrl");

int main(int argc, const char *argv[]) {
  std::cout << "LehrFEM++ Demo program for the use of static control variables"
            << std::endl;
  std::cout << "Variables can be set via command line arguments: VARNAME=VALUE"
            << std::endl;
  std::cout << "Try: testvar=42" << std::endl << std::endl;

  std::cout << ">> List of control variables before initialization:"
            << std::endl;
  lf::base::ListCtrlVars(std::cout);

  // Set verbosity of output for functions called subsequently
  lf::base::SetCtrlVar("read_ctrl_vars_file", 1);
  lf::base::SetCtrlVar("read_ctrl_vars_args", 10);

  lf::base::ReadCtrVarsCmdArgs(argc, argv);

  if (!lf::base::ReadCtrlVarsFile("setup.vars")) {
    std::cout << "No file specifyng control variables" << std::endl;
  }

  std::cout << ">> List of control variables after initialization:"
            << std::endl;
  lf::base::ListCtrlVars(std::cout);

  return 0;
}
