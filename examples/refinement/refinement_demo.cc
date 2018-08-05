/** @file refinement_demo.cc
 * @brief example code for building and refining a mesh
 */

#include <lf/refinement/refinement_hierarchy.h>
#include <lf/refinement/refutils.h>
#include <iostream>
#include "lf/base/base.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

int main() {
  std::cout << "LehrFEM++ demo of mesh construction and refinement"
            << std::endl;

  // Set control variables from command line or file "setup vars"
  // const int n_clvars = lf::base::ReadCtrVarsCmdArgs(argc, argv);
  if (!lf::base::ReadCtrlVarsFile("setup.vars")) {
    std::cout << "No file specifyng control variables" << std::endl;
  }

  return 0;
}
