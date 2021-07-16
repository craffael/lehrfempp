/**
 * @file
 * @brief Visualizes all the (high-order) shape functions of a triangle using
 * high-order vtk files.
 * @author Raffael Casagrande
 * @date   2021-07-10 06:14:24
 * @copyright MIT License
 */

#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

int main() {
  lf::mesh::hybrid2d::MeshFactory mesh_factory(2);
  mesh_factory.AddPoint(Eigen::Vector2d{0, 0});
  mesh_factory.AddPoint(Eigen::Vector2d{3, 0});
  mesh_factory.AddPoint(Eigen::Vector2d{2, 1});
  mesh_factory.AddEntity(lf::base::RefEl::kTria(),
                         std::array<lf::base::size_type, 3>{0, 1, 2}, nullptr);

  auto mesh = mesh_factory.Build();
  const auto* edge0 = mesh->EntityByIndex(1, 0);
  const auto* edge1 = mesh->EntityByIndex(1, 1);
  const auto* edge2 = mesh->EntityByIndex(1, 2);
  const auto* tria = mesh->EntityByIndex(0, 0);

  auto degree_mds = lf::mesh::utils::AllCodimMeshDataSet<unsigned>(mesh);
  degree_mds(*edge0) = 1;
  degree_mds(*edge1) = 2;
  degree_mds(*edge2) = 3;
  degree_mds(*tria) = 3;

  auto fe_space = std::make_shared<lf::fe::HierarchicScalarFESpace<double>>(
      mesh, degree_mds);
  auto coeff_vec = Eigen::VectorXd(fe_space->LocGlobMap().NumDofs());
  coeff_vec.setZero();

  lf::io::VtkWriter vtk(mesh, "out.vtk", 0, 5);

  // Because a ParaView bug
  // (https://gitlab.kitware.com/paraview/paraview/-/issues/20837) we write a
  // complicated function as the first dataset so that the tesselate filter
  // works
  vtk.WritePointData(
      "trig", lf::mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return std::sin(x.x()) * std::sin(x.y());
      }));

  // Write out the individual basis functions
  for (int i = 0; i < coeff_vec.rows(); ++i) {
    coeff_vec(i) = 1.;
    vtk.WritePointData(std::to_string(i),
                       lf::fe::MeshFunctionFE(fe_space, coeff_vec));
    coeff_vec(i) = 0.;
  }
}