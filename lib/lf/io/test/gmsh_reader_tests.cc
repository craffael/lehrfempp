/**
 * @file
 * @brief Test the GmshReader
 * @author Raffael Casagrande
 * @date   2018-07-01 08:11:24
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include "lf/mesh/hybrid2d/hybrid2d.h"

namespace lf::io::test {

using size_type = mesh::Mesh::size_type;

/// Retrieves the name of a mesh that lies in the lib/lf/io/test/msh_files/
/// folder
std::string getMeshPath(std::string mesh_name) {
  std::string file_path = __FILE__;
  // get directory: get last index of / or backslash
  int index = 0;
  char separator;
  for (unsigned int i = 0; i < file_path.length(); ++i) {
    if (file_path[i] == '/' || file_path[i] == '\\') {
      index = i;
      separator = file_path[i];
    }
  }
  std::string directory_path = file_path.substr(0, index + 1);
  return directory_path + "msh_files" + separator + mesh_name;
}

TEST(lf_io, readTwoElementMesh) {
  GmshReader reader(std::make_unique<mesh::hybrid2d::MeshFactory>(2),
                    getMeshPath("two_element_hybrid_2d.msh"));
  auto mesh = reader.mesh();
  EXPECT_EQ(mesh->Size(0), 2);
  EXPECT_EQ(mesh->Size(1), 6);
  EXPECT_EQ(mesh->Size(2), 5);

  auto entities2 = mesh->Entities(2);
  auto origin = std::find_if(entities2.begin(), entities2.end(), [](auto& e) {
    return e.Geometry()->Global(Eigen::MatrixXd(0, 1)).squaredNorm() < 1e-10;
  });
  EXPECT_NE(origin, entities2.end());

  EXPECT_EQ(reader.physicalEntityNr(*origin).size(), 2);
  EXPECT_EQ(reader.physicalEntityNr(*origin)[0], 1);
  EXPECT_EQ(reader.physicalEntityNr(*origin)[1], 2);
}
}  // namespace lf::io::test
