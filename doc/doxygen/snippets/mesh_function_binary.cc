/**
 * @file
 * @brief Shows usage of MeshFunctionBinary and related operator overloads.
 * @author Raffael Casagrande
 * @date   2019-03-02 06:40:56
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe {

void addition() {
  //! [one_trig]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");

  // a mesh function that takes the value 1 everywhere
  auto mf_one = mesh::utils::MeshFunctionConstant(1.);

  // a mesh function which represents the function `sin(x)*cos(y)` (x,y are
  // global coordinates)
  auto mf_trig = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  // mesh function `1+sin(x)*cos(y)`
  auto mf_one_trig = mf_one + mf_trig;
  //! [one_trig]
}

void subtraction() {
  //! [subtract]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");

  // a mesh function that takes the value 1 everywhere
  auto mf_one = mesh::utils::MeshFunctionConstant(1.);

  // a mesh function which represents the function `sin(x)*cos(y)` (x,y are
  // global coordinates)
  auto mf_trig = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  // mesh function `1-sin(x)*cos(y)`
  auto mf_one_minus_trig = mf_one - mf_trig;
  //! [subtract]
}

void multiplication() {
  //! [product]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");

  // a mesh function which takes the value 1 everywhere
  auto mf_one = mesh::utils::MeshFunctionConstant(1.);

  // a mesh function which represents the radial vector field (x,y)
  auto mf_radial = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return x; });

  // a matrix valued mesh function ((x,y),(x^2,y^2))
  auto mf_matrix =
      mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
        return (Eigen::Matrix2d() << x[0], x[1], x[0] * x[0], x[1] * x[1])
            .finished();
      });

  // product of two scalar valued mesh functions:
  auto p0 = mf_one * mesh::utils::MeshFunctionConstant(3.);

  // product of a scalar valued mesh function with a matrix valued one:
  auto p1 = mf_one * mf_matrix;

  // product of matrix valued with a vector valued mesh function:
  auto p2 = mf_matrix * mf_radial;  // will be vector valued

  // product of a matrix valued mesh function with itself:
  auto p3 = mf_matrix * mf_matrix;
  //! [product]
}

}  // namespace lf::uscalfe
