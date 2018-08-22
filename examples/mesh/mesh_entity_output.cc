/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>
#include <vector>
#include <algorithm> // for std::find


#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"
#include "lf/mesh/utils/utils.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/print_info.h"
#include "lf/mesh/utils/utils.h"

int main() {
  
  using namespace lf::mesh::utils;

  std::cout << "Output of information for mesh entity elements" << std::endl;

  // Build mesh ----------------------------------
  lf::mesh::hybrid2dp::MeshFactory test(2);  // MeshFactory object

  // add nodes
  test.AddPoint(Eigen::Vector2d(0, 0));
  test.AddPoint(Eigen::Vector2d(0, 1));
  test.AddPoint(Eigen::Vector2d(1, 1));
  test.AddPoint(Eigen::Vector2d(1, 0));

  // add an element
  Eigen::MatrixXd node_coord(2, 4);
  node_coord << 0, 1, 1, 0, 0, 0, 1, 1;
  test.AddEntity(lf::base::RefEl::kQuad(), {0, 1, 2, 3},
                 std::make_unique<lf::geometry::QuadO1>(std::move(node_coord)));

  // explicitly add the right edge:
  node_coord = Eigen::MatrixXd(2, 2);
  node_coord << 1, 1, 0, 1;
  test.AddEntity(lf::base::RefEl::kSegment(), {1, 2},
                 std::make_unique<lf::geometry::SegmentO1>(node_coord));

  
  // build the mesh and retrieve a pointer
  auto mesh = test.Build();  // mesh is Mesh object // test is MeshFactory object


/*
  // Output information on mesh
  std::cout << "##### Mesh information ######" << std::endl;
  lf::mesh::utils::PrintInfo(*mesh, std::cout);
  std::cout << "#####                   #####" << std::endl;
*/


 // lf::mesh::Entity::output_ctrl_ = 100;


/*
  std::cout << "****** Output of mesh entities *******" << std::endl;
  // Loop over entities and print associated information
  for (lf::base::dim_t codim = 0; codim <= 2; ++codim) {
    std::cout << "******* Entities of codimension " << static_cast<int>(codim)
              << " ******* " << std::endl;
    for (const lf::mesh::Entity &entity : mesh->Entities(codim)) {
      std::cout << entity << std::endl;
    }
  }
*/

  // Test mesh ---------------------------------------------------
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Testing of writeTikZ for mesh (test mesh)
  lf::mesh::utils::writeTikZ(*mesh_p, "tikz_mesh_test2.txt");

/*
                             TikzOutputCtrl::RenderCells|TikzOutputCtrl::CellNumbering|TikzOutputCtrl::VerticeNumbering|
                             TikzOutputCtrl::EdgeNumbering|TikzOutputCtrl::NodeNumbering|TikzOutputCtrl::ArrowTips);
*/




  // write_tikz version 2.0 --------------------------------------

  // List of entities
  std::vector<lf::mesh::Entity> entitiesToPrint = {};

  // Q: How should the list of entities look?
  // Q: How to look through list? std::find?
  // Q: Add global object index?

  auto desiredEntities = [&] (const lf::mesh::Entity& entity){
      //if(std::find(entitiesToPrint.begin(), entitiesToPrint.end(), entity) != entitiesToPrint.end()){
      // if (entitiesToPrint.contains(entity)){return true;} else {return false;}

      if (entitiesToPrint.size() < 2) {
          //std::cout << "true" << std::endl;
          return true;
      } else {
          //std::cout << "false" << std::endl;
          return false;
      }
  }; //auto

  int output_ctrl = 0; // Set default argument
  lf::mesh::utils::writeTikZ2(*mesh_p, "tikz_selector.txt", desiredEntities);

  return 0L;
}
