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
#include <sstream>

#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/refutils.h>
#include "lf/base/base.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

// Global variable for selection of refinement type
// 0 -> regular global refinement
// 1 -> global baryccentric refinement
// 2 -> global marking of edges
// 3 -> local marking of edges with midpoints in [0,1]^2
CONTROLDECLAREINFO(refselector, "refselector",
                   "Selector for refinement method");

int main(int argc, const char *argv[]) {
  using size_type = lf::base::size_type;

  std::cout << "LehrFEM++ demo of mesh construction and refinement"
            << std::endl;

  // Set control variables from command line or file "setup vars"
  lf::base::ReadCtrVarsCmdArgs(argc, argv);
  if (!lf::base::ReadCtrlVarsFile("setup.vars")) {
    std::cout << "No file specifyng control variables" << std::endl;
  }
  std::cout << "##### Control variables:" << std::endl;
  lf::base::ListCtrlVars(std::cout);

  // Generate hybrid test mesh and obtain a pointer to it
  std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  const lf::mesh::Mesh &mesh = *mesh_ptr;

  // Output information about the mesh
  std::cout << "########################################" << std::endl;
  std::cout << "Output of PrintInfo() for test mesh" << std::endl;
  lf::mesh::utils::PrintInfo(mesh, std::cout);

  // Build mesh hierarchy
  std::cout << "########################################" << std::endl;
  std::cout << "Initialization of data structure for mesh hierarchy"
            << std::endl;
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_ptr, mesh_factory_ptr);

  // lambda functions for (local) marking
  // First, a global marker
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &)>
      allmarker = [](const lf::mesh::Mesh &mesh,
                     const lf::mesh::Entity &edge) -> bool { return true; };
  // Mark edges whose center lies inside a square
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &edge)>
      locmarker =
          [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) -> bool {
    Eigen::MatrixXd ref_c(1, 1);
    ref_c(0, 0) = 0.5;
    Eigen::VectorXd c(edge.Geometry()->Global(ref_c));
    return ((c[0] >= 0.0) && (c[0] < 1.0) && (c[1] >= 0.0) && (c[1] < 1.0));
  };

  // Main refinement loop
  std::cout << "########################################" << std::endl;
  std::size_t Nrefs;
  std::cout << "No of refinement steps: ";
  std::cin >> Nrefs;
  std::cout << "Entering main refinement loop" << std::endl;
  for (int refstep = 0; refstep < Nrefs; refstep++) {
    std::cout << "#### Refinement step " << refstep + 1 << std::endl;
    // Depending on the value of the control variable do different types
    // of refinement
    switch (refselector) {
      case 0: {
        // Global regular refinement
        std::cout << "#### global regular refinement" << std::endl;
        multi_mesh.RefineRegular();
        break;
      }
      case 1: {
        // Global barycentric refinement
        std::cout << "#### global barycentric refinement" << std::endl;
        multi_mesh.RefineRegular(lf::refinement::RefPat::rp_barycentric);
        break;
      }
      case 2: {
        // Global regular refinement
        std::cout << "#### global bisection refinement" << std::endl;
        multi_mesh.MarkEdges(allmarker);
        multi_mesh.RefineMarked();
        break;
      }
      case 3: {
        // Global regular refinement
        std::cout << "#### local bisection refinement" << std::endl;
        multi_mesh.MarkEdges(locmarker);
        multi_mesh.RefineMarked();
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Unknown refinement selector");
        break;
      }
    }  // end switch refselector

    // Inquire about number of existing levels
    const size_type n_levels = multi_mesh.NumLevels();
    // Obtain pointer to mesh on finest level
    std::shared_ptr<const lf::mesh::Mesh> mesh =
        multi_mesh.getMesh(n_levels - 1);
    // Print number of entities of various co-dimensions
    std::cout << "#### Mesh on level " << n_levels - 1 << ": " << mesh->Size(2)
              << " nodes, " << mesh->Size(1) << " nodes, " << mesh->Size(0)
              << " cells," << std::endl;
    std::stringstream level_asc;
    level_asc << refstep;
    lf::mesh::utils::writeTikZ(
        *mesh, std::string("refinement_mesh") + level_asc.str() + ".txt", 31);
  }

  // Generate  MATLAB functions that provide a description of all
  // levels of the mesh hierarchy
  std::string basename;
  std::cout << "Basename for MATLAB output: ";
  std::cin >> basename;
  lf::refinement::WriteMatlab(multi_mesh, basename);

  return 0;
}
