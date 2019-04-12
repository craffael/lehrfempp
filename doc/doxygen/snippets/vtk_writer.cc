/**
 * @file
 * @brief Illustrate usage of VtkWriter
 * @author Raffael Casagrande
 * @date   2018-08-13 08:16:01
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::io {

void foo() {
  //! [usage]
  std::shared_ptr<mesh::Mesh> mesh;  // initialize mesh somehow

  // data stored with codim=0 entities
  std::shared_ptr<mesh::utils::MeshDataSet<double>> cell_data;

  // vector data stored with points of the mesh
  std::shared_ptr<mesh::utils::MeshDataSet<Eigen::Vector2d>> point_data;

  // A Mesh function representing the function sin(x)*cos(y)
  auto mf = uscalfe::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  io::VtkWriter vtk_writer(mesh, "filename.vtk");
  vtk_writer.WriteCellData("cellData", *cell_data);
  vtk_writer.WritePointData("pointData", *point_data);
  vtk_writer.WritePointData("mfPoint", mf);
  //! [usage]
}

void mfPointUsage() {
  //! [mfPointUsage]
  std::shared_ptr<mesh::Mesh> mesh;  // initialize mesh somehow

  // construct a first order lagrange fe space
  auto fes = std::make_shared<uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  // A Mesh function representing the first shape function of fes
  Eigen::VectorXd x(fes->LocGlobMap().NoDofs());
  x[0] = 1.0;
  auto mfShapeFun = uscalfe::MeshFunctionFE(fes, x);

  // A Meshfunction representing the gradient of the first shape function:
  auto mfGradShapeFun = uscalfe::MeshFunctionGradFE(fes, x);

  // write to VTK
  io::VtkWriter vtk_writer(mesh, "filename.vtk");
  vtk_writer.WritePointData("shapeFun", mfShapeFun);

  // Trying to sample the gradient directly on the points of the mesh
  // will fail an assert, see note above:
  // vtk_writer.WritePointData("GradShapeFun", mfGradShapeFun);
  //! [mfPointUsage]
}

void mfCellUsage() {
  //! [mfCellUsage]
  std::shared_ptr<mesh::Mesh> mesh;  // initialize mesh somehow

  // construct a first order lagrange fe space
  auto fes = std::make_shared<uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  // A Mesh function representing the gradient of the first shape function
  Eigen::VectorXd x(fes->LocGlobMap().NoDofs());
  x[0] = 1.0;
  auto mfShapeFun = uscalfe::MeshFunctionGradFE(fes, x);

  // write the gradient to vtk (cell based)
  io::VtkWriter vtk_writer(mesh, "filename.vtk");
  vtk_writer.WriteCellData("shapeFun", mfShapeFun);
  //! [mfCellUsage]
}

}  // namespace lf::io
