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

void mfgdemo(const lf::mesh::Mesh& mesh) {
  //! [mfgdemo]
  // Define a tensor field through a lambda function
  auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    return (Eigen::Matrix<double, 2, 2>() << (3.0 + x[1]), x[0], x[0],
            (2.0 + x[0]))
        .finished();
  };
  // Wrap the tensor field into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
  // evaluate the tensor field at the nodes of the mesh
  const Eigen::Matrix<double, 0, 1> dummy;
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    const std::vector<Eigen::Matrix<double, 2, 2>> a{mf_alpha(*node, dummy)};
    std::cout << a[0] << std::endl;
  }
  //! [mfgdemo]
}

//! [mffedemo]
double computeL2ErrorNorm(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    const Eigen::VectorXd& mu) {
  // Reference to underlying mesh
  const lf::mesh::Mesh& mesh{*(fe_space_p->Mesh())};
  // Lambda function providing the "exact solution"
  auto u = [](Eigen::Vector2d x) -> double {
    return std::log(x[0] * x[0] + x[1] + 1.0);
  };
  // Has to be wrapped into a mesh function for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  // create mesh functions representing solution / gradient of solution
  auto mf_sol = lf::uscalfe::MeshFunctionFE(fe_space_p, mu);
  // compute errors with 10-th order quadrature rules
  double L2err_2 =  // NOLINT
      std::sqrt(IntegrateMeshFunction(mesh, squaredNorm(mf_sol - mf_u), 2));
  return std::sqrt(L2err_2);
}
//! [mffedemo]

}  // namespace lf::uscalfe
