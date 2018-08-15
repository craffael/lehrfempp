/**
 * @file
 * @brief Illustrate usage of VtkWriter
 * @author Raffael Casagrande
 * @date   2018-08-13 08:16:01
 * @copyright MIT License
 */

#include <lf/io/io.h>

namespace lf::io {

void foo() {
  //! [usage]
  std::shared_ptr<mesh::Mesh> mesh;  // initialize mesh somehow

  // data stored with codim=0 entities
  std::shared_ptr<mesh::utils::MeshDataSet<double>> cell_data;

  // vector data stored with points of the mesh
  std::shared_ptr<mesh::utils::MeshDataSet<Eigen::Vector2d>> point_data;

  io::VtkWriter vtk_writer(mesh, "filename.vtk");
  vtk_writer.WriteCellData("cellData", *cell_data);
  vtk_writer.WritePointData("pointData", *point_data);
  //! [usage]
}

}  // namespace lf::io
