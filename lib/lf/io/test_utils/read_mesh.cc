/**
 * @file
 * @brief Implementations from read_mesh.h
 * @author Raffael Casagrande
 * @date   2018-07-21 12:01:40
 * @copyright MIT License
 */

#include "read_mesh.h"

namespace lf::io::test_utils {
std::string getMeshPath(std::string mesh_name) {
  std::string file_path = __FILE__;
  // std file_path(__FILE__);
  // get directory: go up two directories last index of / or backslash
  int index = 0;
  char separator;
  char count = 2;
  for (int i = file_path.length() - 1; i >= 0; --i) {
    if (file_path[i] == '/' || file_path[i] == '\\') {
      --count;
      if (count == 0) {
        index = i;
        separator = file_path[i];
        break;
      }
    }
  }
  std::string directory_path = file_path.substr(0, index + 1);
  return directory_path + "test" + separator + "msh_files" + separator +
         mesh_name;
}

GmshReader getGmshReader(std::string mesh_name, base::dim_t dim_world) {
  return GmshReader(std::make_unique<mesh::hybrid2d::MeshFactory>(dim_world),
                    getMeshPath(mesh_name));
}
}  // namespace lf::io::test_utils
