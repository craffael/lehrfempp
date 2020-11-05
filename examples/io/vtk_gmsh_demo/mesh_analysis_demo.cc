/** @file
 * @brief Prints information about a mesh stored in the Gmsh .msh file
 * @author Ralf Hiptmair
 * @date July 2020
 * @copyright MIT License
 */
#include <filesystem>

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <sstream>

int main(int argc, char **argv) {
  using size_type = lf::base::size_type;
  using dim_t = lf::base::dim_t;

  // Expects at most one command line argument, which should be a raw file name
  std::string filename{};
  switch (argc) {
    case 1: {
      // Set default file name
      filename = "bentwire.msh";
      break;
    }
    case 2: {
      // Filename given via the command line
      filename = argv[1];
      break;
    }
    default: {
      std::cerr << "Usage: " << argv[0] << " <filename> " << std::endl;
      break;
    }
  }  // end switch

  // Mesh file is supposed to reside in the same directory as source file
  std::filesystem::path here = __FILE__;
  auto mesh_path = here.parent_path() / filename.c_str();

  // load the mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());
  std::cout << "Mesh read from file " << mesh_path.string() << std::endl;

  // Obtain pointer to mesh object
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{reader.mesh()};
  const lf::mesh::Mesh &mesh{*mesh_p};
  std::cout << mesh.DimMesh() << "D mesh in " << mesh.DimWorld()
            << "D: " << mesh.NumEntities(2) << " nodes, " << mesh.NumEntities(1)
            << " edges, " << mesh.NumEntities(0) << " cells" << std::endl;

  // Examine physical groups
  for (dim_t codim = 0; codim <= 2; codim++) {
    // Obtain list of physical groups for that co-dimension
    std::vector<std::pair<size_type, std::string>> phys_ent_list{
        reader.PhysicalEntities(codim)};
    if (!phys_ent_list.empty()) {
      // Loop over all physical groups of a particular co-dimension
      for (auto &ent_code : phys_ent_list) {
        std::cout << "codim = " << codim << " physical group "
                  << ent_code.second << ", id = " << ent_code.first;
        // Count mesh entities belonging to that physical group
        int cnt{0};
        for (const lf::mesh::Entity *entity : mesh.Entities(codim)) {
          if (reader.IsPhysicalEntity(*entity, ent_code.first)) {
            cnt++;
          }
        }
        std::cout << ": " << cnt << " entitites" << std::endl;
      }
    } else {
      std::cout << "codim = " << codim << ": No physical groups" << std::endl;
    }
  }

  // Finally output (small) mesh in TikZ format: basic information
  {
    using lf::io::TikzOutputCtrl;
    std::stringstream filename_tikz;
    filename_tikz << filename.substr(0, filename.find_last_of('.')) << ".tex";
    lf::io::writeTikZ(*mesh_p, filename_tikz.str(),
                      TikzOutputCtrl::WithPreamble);
  }
  // Finally output (small) mesh in TikZ format: full information
  {
    using lf::io::TikzOutputCtrl;
    std::stringstream filename_tikz;
    filename_tikz << filename.substr(0, filename.find_last_of('.'))
                  << "_full.tex";
    lf::io::writeTikZ(
        *mesh_p, filename_tikz.str(),
        TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
            TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
            TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::ArrowTips |
            TikzOutputCtrl::WithPreamble);
  }
}  // end main

/* Output when invoked without an argument
tk_gmsh_demo $ ./examples.io.mesh_analysis_demo
Mesh read from file
/scratch/users/ralfh/NOSAVE/numpde/Numcourses/NumPDE/lehrfempp/examples/io/vtk_gmsh_demo/bentwire.msh
2D mesh in 2D: 55 nodes, 131 edges, 77 cells
codim = 0 physical group wire, id = 4: 77 entitites
codim = 1 physical group Contact0, id = 1: 2 entitites
codim = 1 physical group Contact1, id = 2: 2 entitites
codim = 1 physical group Insulated, id = 3: 27 entitites
codim = 2 physical group bottomcontact, id = 5: 2 entitites
codim = 2 physical group leftcontact, id = 6: 2 entitites
*/
