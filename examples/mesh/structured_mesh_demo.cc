/** @brief Demonstrates the construction of a structured mesh in LehrFEM++ */

#include <boost/program_options.hpp>
#include <iostream>

#include "lf/base/base.h"
#include "lf/io/io.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

int main(int argc, char** argv) {
  using size_type = lf::base::size_type;

  namespace po = boost::program_options;

  // clang-format off
  po::options_description desc("LehrFEM++ demo: construction of tensor producy triangular mesh\nUse the option `-h` to display the control variables.");
  desc.add_options()
    ("help,h", "Display help")
    ("Nx_cells", po::value<int>()->required(), "No of cells in X direction, must be > 0")
    ("Ny_cells", po::value<int>()->required(), "No of cells in Y direction, must be > 0")
  ;
  //clang-format on

  // Set control variables from command line or file "setup vars"
  // check for file with options
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  try {
    po::store(po::parse_config_file<char>("setup.vars", desc), vm);
  } catch(const po::reading_file& error) {
    std::cout << "No file `setup.vars` specifying control variables\n";
  }

  if(vm.count("help") > 0) {
    std::cout << desc << '\n';
    return 1;
  }

  po::notify(vm);

  

  // number of cells cannot be 0, that would result in segfault!
  const int Nx = vm["Nx_cells"].as<int>();
  const int Ny = vm["Ny_cells"].as<int>();
  if (Nx == 0 || Ny == 0) {
    std::cout << "Nx and Ny must not be zero, set using options --Nx_cells and "
                 "--Ny_cells\n";
    return 0;
  }

  // Construct a triangular tensor product mesh with 2*Nx*Ny cells
  lf::mesh::utils::TPTriagMeshBuilder builder(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(Nx)
      .setNumYCells(Ny);
  auto mesh_p = builder.Build();

  // Some consistency checks.
  // These should not be included in a regular code
  std::cout << "Checking entity indexing" << '\n';
  lf::mesh::test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << '\n';
  lf::mesh::test_utils::checkMeshCompleteness(*mesh_p);

  // Printing mesh information
  lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

  // Matlab output of mesh
  lf::io::writeMatlab(*mesh_p, "lf_tptriagmesh.m");
  return 0;
}
