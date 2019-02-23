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

ADDOPTION(Nx, Nx_cells, "No of cells in X direction, must be > 0");
ADDOPTION(Ny, Ny_cells, "No of cells in Y direction, must be > 0");

int main(int argc, char** argv) {
  using size_type = lf::base::size_type;

  std::cout << "LehrFEM++ demo: construction of tensor producy triangular mesh\n"
            << "Use the option `-h` to display the control variables.\n";

  // Set control variables from command line or file "setup vars"
  lf::base::ci::Add("help,h", "Display help");
  // check for file with options 
  if (!lf::base::ci::ParseFile("setup.vars")) {
    std::cout << "No file `setup.vars` specifying control variables\n";
  }
  lf::base::ci::ParseCommandLine(argc, argv);
  if (lf::base::ci::Help())
    return 0;

  // number of cells cannot be 0, that would result in segfault!
  if (Nx == 0 || Ny == 0) {
    std::cout << "Nx and Ny must not be zero, set using options --Nx_cells and --Ny_cells\n";
    return 0;
  }

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
