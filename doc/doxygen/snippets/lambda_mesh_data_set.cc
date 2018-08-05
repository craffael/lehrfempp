/**
 * @file
 * @brief Shows usage of LambdaMeshDataSet.
 * @author Raffael Casagrande
 * @date   2018-08-04 04:51:52
 * @copyright MIT License
 */

#include "lf/mesh/utils/lambda_mesh_data_set.h"
#include <lf/io/test_utils/read_mesh.h>

namespace lf::mesh::utils {

void foo() {
  //! [usage]
  // MeshDataSet that stores a `1` with every entity of the mesh
  auto ones_mds = make_LambdaMeshDataSet([](const Entity& e) { return 1; });

  // A MeshDataSet that stores the codimension of every entity:
  auto index_mds =
      make_LambdaMeshDataSet([](const auto& e) { return e.Codim(); });

  // A MeshDataSet that stores the string "hello" with every codim=0 entity
  // and that is undefined for other codimensions:
  auto hello_mds =
      make_LambdaMeshDataSet([](const auto& e) { return "hello"; },
                             [](const auto& e) { return e.Codim() == 0; });

  //! [usage]
}

}  // namespace lf::mesh::utils
