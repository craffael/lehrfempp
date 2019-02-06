/** @file
 * @brief Output of test meshes
 * @author Ralf Hiptmair
 * @date December 2018
 * @copyright MIT License
 */

#include <sstream>
#include "lf/io/io.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

int main() {
  using lf::io::TikzOutputCtrl;
  std::cout << "LehrFEM++ demo: output of test meshes" << std::endl;
  std::cout << "(test meshes from "
               "lf::mesh::test_utils::GenerateHybrid2DTestMesh()"
            << std::endl;
  for (int selector = 0;
       selector <= lf::mesh::test_utils::GenerateHybrid2DTestMesh_maxsel;
       selector++) {
    // Generates a small test meshes, more precisely described in the
    // documentation of the function GenerateHybrid2DTestMesh()
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    // Output of mesh information
    std::cout << "#### Test mesh " << selector << " ####" << std::endl;
    std::cout << *mesh_p << std::endl;

    // TikZ output
    std::stringstream filename_tikz;
    filename_tikz << "test_mesh_" << selector << ".tex";
    lf::io::writeTikZ(
        *mesh_p, filename_tikz.str(),
        TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
            TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
            TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::ArrowTips);

    // Matlab output
    std::stringstream filename_matlab;
    filename_matlab << "test_mesh_" << selector << ".m";
    lf::io::writeMatlab(*mesh_p, filename_matlab.str());

    // Python output
    std::stringstream filename_py;
    filename_py << "test_mesh_" << selector << ".py";
    lf::io::writeMatplotlib(*mesh_p, filename_py.str());

    // VTK output
    std::stringstream filename_vtk;
    filename_vtk << "test_mesh_" << selector << ".vtk";
    lf::io::VtkWriter vtk_writer(mesh_p, filename_vtk.str());
  }

  return 0L;
}
