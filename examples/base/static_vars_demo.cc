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
  static unsigned int ctrl_var_;
  static unsigned int other_var_;
  static unsigned int arg_var_;
  static unsigned int output_ctrl_;  // Added
};

// Declaration of static control variables. They are all of type int.
CONTROLDECLARE(testvar, "testvar");
// Declare StaticVarsDemoClass::ctrl_var_ and set it to 0.
CLASSCONTROLDECLARE(StaticVarsDemoClass, ctrl_var_, "ctrl_var");
// Same as above, but add a comment to the variable. We can later print
// this comment to inform the user about the variable.
CONTROLDECLARECOMMENT(StaticVarsDemoClass, other_var_, "other_var",
                      "A test case");

CONTROLDECLARECOMMENT(StaticVarsDemoClass, arg_var_, "arg_var",
                      "Set from command line");

CONTROLDECLARECOMMENT(StaticVarsDemoClass, output_ctrl_, "output_ctrl_",
                      "Test output ctrl");

int main(int argc, const char *argv[]) {
  std::cout << "LehrFEM++ Demo program for the use of static control variables."
            << "\n"
            << "Variables can be set via command line arguments: VARNAME=VALUE"
            << " (but they must be declared in the program).\n"
            << "Try: testvar=42 \n\n";

  std::cout << ">> List of control variables before initialization:"
            << "\n";
  lf::base::ListCtrlVars(std::cout);
  std::cout << "\n"; // for sake of readability

  // Set verbosity of output for functions called subsequently
  // The static variables `read_ctrl_vars_file` and `read_ctrl_vars_args`
  // are already declared in static_vars.cc.
  // `read_ctrl_vars_file`: if > 0, reading variables from files is enabled
  // `read_ctrl_vars_args`: if > 0, command line arguments will be parsed for
  //                        static variables
  lf::base::SetCtrlVar("read_ctrl_vars_file", 1);
  lf::base::SetCtrlVar("read_ctrl_vars_args", 1); // 10? couln't it just be 1?

  // Parse argv for definitions of the static variables we declared above
  lf::base::ReadCtrVarsCmdArgs(argc, argv);

  // Check the file "setup.vars" for the definition of static variables
  if (!lf::base::ReadCtrlVarsFile("setup.vars")) {
    std::cout << "No file specifying control variables.\n\n";
  }

  std::cout << "\n>> List of control variables after initialization:\n";
  lf::base::ListCtrlVars(std::cout);

  return 0;
}
