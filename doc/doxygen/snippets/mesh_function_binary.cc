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
  auto mf_one = MeshFunctionConstant(1.);

  // a mesh function which represents the function `sin(x)*cos(y)` (x,y are
  // global coordinates)
  auto mf_trig = MeshFunctionGlobal(
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
  auto mf_one = MeshFunctionConstant(1.);

  // a mesh function which represents the function `sin(x)*cos(y)` (x,y are
  // global coordinates)
  auto mf_trig = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  // mesh function `1-sin(x)*cos(y)`
  auto mf_one_minus_trig = mf_one - mf_trig;
  //! [subtract]
}

}  // namespace lf::uscalfe
