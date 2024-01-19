/**
 * @file
 * @brief Reads meshes containing second order elements and performs refinement
 * @author Anian Ruoss
 * @date   24.02.2019 18:24:17
 * @copyright MIT License
 */

#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/refinement/mesh_hierarchy.h>

#include <filesystem>
#include <iostream>

using lf::io::TikzOutputCtrl;

int main() {
  const std::filesystem::path file_path = __FILE__;

  for (const std::string& mesh_name :
       {"square_quads.msh", "square_trias.msh"}) {
    auto mesh_path = file_path.parent_path() / "meshes" / mesh_name;

    // read mesh from file
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory),
                                    mesh_path.string());

    // create mesh hierarchy from mesh for refinement
    lf::refinement::MeshHierarchy multi_mesh(
        std::const_pointer_cast<lf::mesh::Mesh>(reader.mesh()),
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));

    const int refinement_steps = 2;

    for (int step = 0; step < refinement_steps; ++step) {
      // refine mesh and store to TikZ
      multi_mesh.RefineRegular();
      auto mesh = multi_mesh.getMesh(multi_mesh.NumLevels() - 1);
      lf::io::writeTikZ(
          *mesh,
          mesh_name.substr(0, mesh_name.find_last_of('.')) + "_" +
              std::to_string(step) + ".tex",
          TikzOutputCtrl::RenderCells | TikzOutputCtrl::CellNumbering |
              TikzOutputCtrl::VerticeNumbering | TikzOutputCtrl::NodeNumbering |
              TikzOutputCtrl::EdgeNumbering);
    }
  }

  return 0;
}
