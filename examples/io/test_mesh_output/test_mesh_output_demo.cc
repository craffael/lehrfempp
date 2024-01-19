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

/**
 * @brief "graphical" output of test meshes in various formats
 * @param mesh_p pointer to LehrFEM++ mesh object to be rendered
 * @param filename name of the file to be written. An appropriate
 * suffix will be added.
 *
 * This function writes the mesh geometry and topology information
 * to file in various formats for graphical rendering:
 * - LaTeX TikZ format, stored in `.tex` file
 * - Matlab format, to be rendered by Matlab function  `plot_lf_mesh.m`
 * - Python format meant for Python script 'plot_mesh.py'
 * - vtk format for rendering with Paraview
 *
 * @note We have to pass pointer to the mesh, because vtk output requires the
 * the construction of a writer object that holds a shared pointer to the mesh.
 */
void writeMeshRenderingData(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
                            const char* filename);

void writeMeshRenderingData(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
                            const char* filename) {
  // TikZ output
  using lf::io::TikzOutputCtrl;
  std::stringstream filename_tikz;
  filename_tikz << filename << ".tex";
  lf::io::writeTikZ(
      *mesh_p, filename_tikz.str(),
      TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
          TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
          TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::ArrowTips);

  // Matlab output
  std::stringstream filename_matlab;
  filename_matlab << filename << ".m";
  lf::io::writeMatlab(*mesh_p, filename_matlab.str());

  // Python output
  std::stringstream filename_py;
  filename_py << filename;
  lf::io::writeMatplotlib(*mesh_p, filename_py.str());

  // VTK output
  std::stringstream filename_vtk;
  filename_vtk << filename << ".vtk";
  const lf::io::VtkWriter vtk_writer(mesh_p, filename_vtk.str());
}

int main() {
  std::cout << "LehrFEM++ demo: output of test meshes" << '\n';
  std::cout << "(test meshes from "
               "lf::mesh::test_utils::GenerateHybrid2DTestMesh()"
            << '\n';
  for (int selector = 0;
       selector <= lf::mesh::test_utils::GenerateHybrid2DTestMesh_maxsel;
       selector++) {
    // Generates a small test meshes, more precisely described in the
    // documentation of the function GenerateHybrid2DTestMesh()
    const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    // Output of mesh information
    std::cout << "#### Test mesh " << selector << " ####" << '\n';
    lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

    std::stringstream filename;
    filename << "test_mesh_" << selector;
    writeMeshRenderingData(mesh_p, filename.str().c_str());
  }

  {
    // Initialize builder object
    lf::mesh::utils::TPTriagMeshBuilder builder(
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
    // Set mesh parameters following the Builder pattern
    // Domain is the unit square, two cells in every direction
    builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
        .setTopRightCorner(Eigen::Vector2d{1, 1})
        .setNumXCells(2)
        .setNumYCells(2);
    const std::shared_ptr<const lf::mesh::Mesh> mesh_p = builder.Build();
    // Output of mesh information
    std::cout << "#### Triangular tensor product 2x2 mesh #####" << '\n';
    lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

    writeMeshRenderingData(mesh_p, "tp_tria_mesh");
  }
  return 0L;
}
