
#include <cassert>
#include <iostream>

#include "lf/mesh/mesh.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

namespace lf::mesh::utils {
void foo() {
  //![writeTikzUsage]
  // Given a mesh object named mesh
  // With the use of the enum flag for node numbering
  // writeTikZ(*mesh, "filename.txt", TikzOutputCtrl::NodeNumbering);

  // Combining enum flags, enabling more detailed output
  // writeTikZ(*mesh, "filename.txt",
  // TikzOutputCtrl::NodeNumbering|TikzOutputCtrl::EdgeNumbering|TikzOutputCtrl::VerticeNumbering);

  // Without flags

  //![writeTikzUsage]

  //![TikzInLatex]

  // \documentclass{article}
  // \usepackage{tikz}
  // \begin{document}

  // \input{"filename.txt"}

  // \end{document}
  //![TikzInLatex]

}  // foo_utils

void AllCodimMeshDataSet_Example() {
  //![AllCodimMeshDataSet]
  // Create a mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Create a MeshDataSet that stores a boolean for every entity
  auto mesh_data_set = lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p);

  // Create a MeshDataSet that stores a boolean for every entity with default
  // value
  auto mds_default = lf::mesh::utils::AllCodimMeshDataSet<bool>(mesh_p, true);

  // set the value associated with an entity to false
  auto entity = mesh_p->EntityByIndex(0, 0);
  mesh_data_set(*entity) = false;

  // check if the data set is defined on an entity
  mesh_data_set.DefinedOn(*entity);  // returns true
  //![AllCodimMeshDataSet]
}

void CodimMeshDataSet_Example() {
  //![CodimMeshDataSet]
  // Create a mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Create a MeshDataSet that stores a boolean for every cell
  auto mesh_data_set_cells = lf::mesh::utils::CodimMeshDataSet<bool>(mesh_p, 0);

  // Create a MeshDataSet that stores a boolean for every edge with default
  // value
  auto mds_edge = lf::mesh::utils::CodimMeshDataSet<bool>(mesh_p, 1, true);

  // set the value associated with a cell to false
  auto cell = mesh_p->EntityByIndex(0, 0);
  mesh_data_set_cells(*cell) = false;

  // check if the data set is defined on an entity
  mesh_data_set_cells.DefinedOn(*cell);  // returns true
  //![CodimMeshDataSet]
}

}  // namespace lf::mesh::utils
