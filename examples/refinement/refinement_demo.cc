/** @file refinement_demo.cc
 * @brief example code for building and refining a mesh
 * 
 * This code generates a test mesh and performs a user specified number of 
 * steps of regular uniform refinement. 
 * 
 * In the end the information about all meshes created in the process 
 * is stored in the form of MATLAB functions.
 */

#include <iostream>

#include <lf/refinement/refinement_hierarchy.h>
#include <lf/refinement/refutils.h>
#include <iostream>
#include "lf/base/base.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"

int main() {
  using size_type = lf::base::size_type;
  
  std::cout << "LehrFEM++ demo of mesh construction and refinement"
            << std::endl;

  // Set control variables from command line or file "setup vars"
  // const int n_clvars = lf::base::ReadCtrVarsCmdArgs(argc, argv);
  if (!lf::base::ReadCtrlVarsFile("setup.vars")) {
    std::cout << "No file specifyng control variables" << std::endl;
  }

  // Generate hybrid test mesh and obtain a pointer to it
  std::shared_ptr<lf::mesh::Mesh>
    mesh_ptr = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  const lf::mesh::Mesh &mesh = *mesh_ptr;

  // Output information about the mesh
  std::cout << "########################################" << std::endl;
  std::cout << "Output of PrintInfo() for test mesh" << std::endl;
  lf::mesh::utils::PrintInfo(mesh,std::cout);

  // Build mesh hierarchy
  std::cout << "########################################" << std::endl;
  std::cout << "Initialization of data structure for mesh hierarchy" << std::endl;
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_ptr, mesh_factory_ptr);

  // Main refinement loop
  std::cout << "########################################" << std::endl;
  std::size_t Nrefs; std::cout << "No of refinement steps: "; std::cin >> Nrefs;
  std::cout << "Entering main refinement loop" << std::endl;
  for (int refstep = 0; refstep < Nrefs; refstep++) {
    std::cout << "#### Refinement step " << refstep+1 << std::endl;
    // Global regular refinement
    multi_mesh.RefineRegular();
    // Inquire about number of existing levels
    const size_type n_levels = multi_mesh.NumLevels();
    // Obtain pointer to mesh on finest level
    std::shared_ptr<const lf::mesh::Mesh> mesh = multi_mesh.getMesh(n_levels - 1);
    // Print number of entities of various co-dimensions
    std::cout << "#### Mesh on level " << n_levels - 1 << ": " << mesh->Size(2)
              << " nodes, " << mesh->Size(1) << " nodes, " << mesh->Size(0)
              << " cells," << std::endl;
  }

  // Generate  MATLAB functions that provide a description of all
  // levels of the mesh hierarchy
  std::string basename; std::cout << "Basename for MATLAB output: ";
  std::cin >> basename;
  lf::refinement::WriteMatlab(multi_mesh,basename);
  
  return 0;
}
