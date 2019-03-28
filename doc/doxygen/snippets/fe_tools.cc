/**
 * @file
 * @brief Doxygens snippets for functions in fe_tools.h
 * @author Raffael Casagrande
 * @date   2019-03-04 10:51:20
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe {

void integrateMeshFunction() {
  //! [integrateMeshFunction]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  // integrate the function sin(x)*cos(y) over the mesh using 5th-degree
  // quadrature rules
  auto mf = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  auto integral = IntegrateMeshFunction(*mesh, mf, 5);
  //! [integrateMeshFunction]
}

void integrateMeshFunction2() {
  //! [integrateMeshFunction2]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  // integrate the function sin(x)*cos(y) over the mesh using 5th-degree
  // quadrature rules
  auto mf = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  // select the quadrature rule explicitly for every element:
  auto integral = IntegrateMeshFunction(*mesh, mf, [](const mesh::Entity& e) {
    return quad::make_QuadRule(e.RefEl(), 5);
  });
  //! [integrateMeshFunction2]
}

void nodalProjection() {
  //! [nodalProjection]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  // project a first-degree polynomial onto a first order lagrange space and
  // make sure the representation is exact:
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);
  auto mf_linear = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return 2 + 3 * x[0] + 4 * x[1]; });

  auto dof_vector = NodalProjection(*fe_space, mf_linear);
  auto mf_fe = MeshFunctionFE(fe_space, dof_vector);

  assert(IntegrateMeshFunction(*mesh, squaredNorm(mf_fe - mf_linear), 2) <
         1e-12);
  //! [nodalProjection]
}

}  // namespace lf::uscalfe
