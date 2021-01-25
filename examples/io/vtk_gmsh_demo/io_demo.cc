/**
 * @file
 * @brief Shows how the classes in the namespace lf::io can be used.
 * @author Raffael Casagrande
 * @date   2018-09-02 02:44:15
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>

#include <filesystem>

int main() {
  // Find path to the smiley mesh
  std::filesystem::path here = __FILE__;
  auto smiley_path = here.parent_path() / "smiley.msh";

  // load the smiley mesh
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), smiley_path.string());

  // print all physical entities:
  std::cout << "Physical Entities in Gmsh File " << std::endl;
  std::cout
      << "---------------------------------------------------------------\n";
  for (lf::base::dim_t codim = 0; codim <= 2; ++codim) {
    for (auto& pair : reader.PhysicalEntities(codim)) {
      std::cout << "codim = " << static_cast<int>(codim) << ": " << pair.first
                << " <=> " << pair.second << std::endl;
    }
  }
  std::cout << std::endl << std::endl;

  // count the number of elements in the eyes:
  int num_eye_elements = 0;
  auto mesh = reader.mesh();
  auto physical_entity_nr_eyes = reader.PhysicalEntityName2Nr("eyes");
  for (auto e : mesh->Entities(0)) {
    if (reader.IsPhysicalEntity(*e, physical_entity_nr_eyes)) {
      ++num_eye_elements;
    }
  }
  std::cout << "Number of eye elements : " << num_eye_elements << std::endl;

  // count the number of edges that form the mouth
  int num_mouth_edges = 0;
  auto physical_entity_nr_mouth = reader.PhysicalEntityName2Nr("mouth");
  for (auto e : mesh->Entities(1)) {
    if (reader.IsPhysicalEntity(*e, physical_entity_nr_mouth)) {
      ++num_mouth_edges;
    }
  }
  std::cout << "Number of mouth edges : " << num_mouth_edges << std::endl;

  // count the total number of elements in the face (including eyes)
  int num_face_elements = 0;
  auto physical_entity_nr_face = reader.PhysicalEntityName2Nr("face");
  for (auto e : mesh->Entities(0)) {
    if (reader.IsPhysicalEntity(*e, physical_entity_nr_face)) {
      ++num_face_elements;
    }
  }
  std::cout << "Number of face elements : " << num_face_elements << std::endl;

  // Start vtk output:
  lf::io::VtkWriter vtk_writer(mesh, "smiley.vtk");

  // write a cell-based function that takes the value 1 on the eyes and 0
  // everywhere else
  {
    lf::mesh::utils::CodimMeshDataSet<int> data(mesh, 0);
    for (auto e : mesh->Entities(0)) {
      data(*e) = static_cast<int>(
          reader.IsPhysicalEntity(*e, physical_entity_nr_eyes));
    }
    vtk_writer.WriteCellData("eyes", data);
  }

  // write nodal data that takes the value 2 on the eyes and 1 on the mouth,
  // zero everywhere else.
  auto mds = lf::mesh::utils::make_CodimMeshDataSet<int>(mesh, 2, 0);
  // mark the eyes
  for (auto e : mesh->Entities(0)) {
    if (reader.IsPhysicalEntity(*e, physical_entity_nr_eyes)) {
      for (auto node : e->SubEntities(2)) {
        mds->operator()(*node) = 2;
      }
    }
  }
  // mark the mouth
  for (auto e : mesh->Entities(1)) {
    if (reader.IsPhysicalEntity(*e, physical_entity_nr_mouth)) {
      for (auto node : e->SubEntities(1)) {
        mds->operator()(*node) = 1;
      }
    }
  }
  vtk_writer.WritePointData("smile", *mds);

  // write a node based, vector valued function that points to the nose:
  vtk_writer.WritePointData(
      "nose_vec", lf::mesh::utils::MeshFunctionGlobal(
                      [&](const Eigen::Vector2d& x) { return (-x).eval(); }));
}
