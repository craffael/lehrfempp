/** @brief Demonstrates the construction of a structured mesh in LehrFEM++ */

#include <iostream>

#include "lf/base/base.h"
#include "lf/io/io.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

CONTROLDECLAREINFO(Nx, "Nx", "No of cells in X direction");
CONTROLDECLAREINFO(Ny, "Ny", "No of cells in X direction");

int main(int argc, const char *argv[]) {
  using size_type = lf::base::size_type;

  std::cout << "LehrFEM++ demo: construction of tensor producy triangular mesh"
            << std::endl;

  // Set control variables from command line or file "setup vars"
  lf::base::ReadCtrVarsCmdArgs(argc, argv);
  if (!lf::base::ReadCtrlVarsFile("setup.vars")) {
    std::cout << "No file specifyng control variables" << std::endl;
  }
  std::cout << "##### Control variables:" << std::endl;
  lf::base::ListCtrlVars(std::cout);

  // Construct a triangular tensor product mesh with 2*Nx*Ny cells
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNoXCells(Nx)
      .setNoYCells(Ny);
  auto mesh_p = builder.Build();

  // Some consistency checks.
  // These should not be included in a regular code
  std::cout << "Checking entity indexing" << std::endl;
  lf::mesh::test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  lf::mesh::test_utils::checkMeshCompleteness(*mesh_p);

  // Printing mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);

  // Matlab output of mesh
  lf::io::writeMatlab(*mesh_p, "lf_tptriagmesh.m");
  return 0;
}
