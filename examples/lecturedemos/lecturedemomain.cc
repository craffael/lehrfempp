/**
 * @file
 * @brief Driver function for simple LehrFEM++ demo
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include "lecturedemomesh.h"
#include <boost/filesystem.hpp>

int main() {
  // Find path to the sample mesh in Gmsh format
  boost::filesystem::path here = __FILE__;
  auto mesh_file = here.parent_path() / "lecturedemomesh.msh";

  // Create a 2D mesh data structure from the information contained in the file
  // A factory object is in charge of creating mesh entities
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file.string());
  // Obtain read only mesh: a view
  std::shared_ptr<const lf::mesh::Mesh> mesh_ptr = reader.mesh();
  const lf::mesh::Mesh &mesh{*mesh_ptr};

  // Output general information on mesh
  std::cout << "Mesh from file " << mesh_file.string() << ": ["
            << (int)mesh.DimMesh() << ',' << (int)mesh.DimWorld()
            << "] dim:" << std::endl;
  std::cout << mesh.NumEntities(0) << " cells, " << mesh.NumEntities(1)
            << " edges, " << mesh.NumEntities(2) << " nodes" << std::endl;

  // First demo: container functionalty of a mesh object
  for (lf::base::dim_t codim = 0; codim <= 2; ++codim) {
    (void)lecturedemo::traverseEntities(mesh, codim);
  }
  // Output information about topological relationships
  lecturedemo::scanTopology(mesh, 0);  // topolgy from the cell perspective
  lecturedemo::scanTopology(mesh, 1);  // topolgy from the edge perspective

  // Output rendering of mesh in LaTeX/TikZ format
  using lf::io::TikzOutputCtrl;
  lf::io::writeTikZ(
      mesh, "demomesh.tex",
      TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
          TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
          TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::ArrowTips);

  return 0L;
}
