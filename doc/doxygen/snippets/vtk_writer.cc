/**
 * @file
 * @brief Illustrate usage of VtkWriter
 * @author Raffael Casagrande
 * @date   2018-08-13 08:16:01
 * @copyright MIT License
 */

#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
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
  auto mf = mesh::utils::MeshFunctionGlobal(
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
  Eigen::VectorXd x(fes->LocGlobMap().NumDofs());
  x[0] = 1.0;
  auto mfShapeFun = fe::MeshFunctionFE(fes, x);

  // A Meshfunction representing the gradient of the first shape function:
  auto mfGradShapeFun = fe::MeshFunctionGradFE(fes, x);

  // write to VTK
  io::VtkWriter vtk_writer(mesh, "filename.vtk");
  vtk_writer.WritePointData("shapeFun", mfShapeFun);

  // Trying to sample the gradient directly on the points of the mesh
  // will fail an assert, see note above:
  // vtk_writer.WritePointData("GradShapeFun", mfGradShapeFun);

  // Sample the geometry and the mesh function with 3rd order lagrange basis
  // functions (see documentation of VtkWriter):
  io::VtkWriter vtk_writerOrder3(mesh, "order3.vtk", 0, 3);
  auto mfTrig = mesh::utils::MeshFunctionGlobal(
      [](const auto& x) { return std::sin(x[0]) * std::cos(x[1]); });
  vtk_writer.WritePointData("mfTrig", mfTrig);

  //! [mfPointUsage]
}

void mfCellUsage() {
  //! [mfCellUsage]
  std::shared_ptr<mesh::Mesh> mesh;  // initialize mesh somehow

  // construct a first order lagrange fe space
  auto fes = std::make_shared<uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  // A Mesh function representing the gradient of the first shape function
  Eigen::VectorXd x(fes->LocGlobMap().NumDofs());
  x[0] = 1.0;
  auto mfShapeFun = fe::MeshFunctionGradFE(fes, x);

  // write the gradient to vtk (cell based)
  io::VtkWriter vtk_writer(mesh, "filename.vtk");
  vtk_writer.WriteCellData("shapeFun", mfShapeFun);
  //! [mfCellUsage]
}

void highOrder() {
  //! [highOrder]
  // read a 2nd order gmsh mesh:
  GmshReader reader(std::make_unique<mesh::hybrid2d::MeshFactory>(2),
                    "unit_circle.msh");
  std::shared_ptr<mesh::Mesh> mesh = reader.mesh();

  // define mesh function of the form sin(\pi*x)*cos(\pi*y) (in global
  // coordinates)
  auto mfTrig = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return std::sin(base::kPi * x[0]) * std::cos(base::kPi * x[1]);
  });

  // output the mesh + mesh function using 1st order cells:
  VtkWriter vtk1(mesh, "1storder.vtk");
  vtk1.WritePointData("trig", mfTrig);

  // output the mesh + mesh function using 2nd order cells:
  VtkWriter vtk2(mesh, "2ndorder.vtk", 0, 2);
  vtk2.WritePointData("trig", mfTrig);
  //! [highOrder]
}

}  // namespace lf::io
