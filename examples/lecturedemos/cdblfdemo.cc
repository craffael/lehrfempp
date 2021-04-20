/**
 * @file
 * @brief demonstration of assembly of Galerkin linear system in LehrFEM++
 * assemble module; meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   April 2021
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

namespace cblfdemo {
// Auxiliary function computing the gradients of barycentric coordinate
// functions on a striangle and returning them in the columns of 2x3 matrix
std::pair<Eigen::Matrix<double, 2, 3>, double> getGradBaryCoords(
    const lf::mesh::Entity& tria) {
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "getGradBeyCoords: Unsupported cell type " << tria.RefEl());
  // Obtain vertex coordinates of the triangle in a 2x3 matrix
  const auto vertices = lf::geometry::Corners(*(tria.Geometry()));
  // Set up an auxiliary 3x3-matrix with a leading column 1 and
  // the vertex coordinates in its right 3x2 block
  Eigen::Matrix<double, 3, 3> X;  // temporary matrix
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = vertices.transpose();
  // The determinant of the auxiliary matrix also supplies the determinant
  // Compute the gradients of the barycentric coordinate functions
  // and store them in the columns of a 2x3 matrix
  return {X.inverse().block<2, 3>(1, 0), 0.5 * std::abs(X.determinant())};
}

/** @brief Compute element matrix for convection-diffusion bilinear form on
 * **triangles**
 */
template <typename MeshFunction>
class CDBLFElemMatProvider {
 public:
  // No default constructor
  CDBLFElemMatProvider() = delete;
  // Constructor takes a vectorfield argument
  explicit CDBLFElemMatProvider(MeshFunction av) : av_(std::move(av)) {}
  // Other constructors
  CDBLFElemMatProvider(const CDBLFElemMatProvider&) = default;
  CDBLFElemMatProvider(CDBLFElemMatProvider&&) noexcept = default;
  CDBLFElemMatProvider& operator=(const CDBLFElemMatProvider&) = delete;
  CDBLFElemMatProvider& operator=(CDBLFElemMatProvider&&) = delete;
  // The crucial interface methods for ENTITY_MATRIX_PROVIDER
  virtual bool isActive(const lf::mesh::Entity& /*cell*/) { return true; }
  Eigen::Matrix<double, 1, 3> Eval(const lf::mesh::Entity& tria);

 private:
  MeshFunction av_;
};

template <typename MeshFunction>
Eigen::Matrix<double, 1, 3> CDBLFElemMatProvider<MeshFunction>::Eval(
    const lf::mesh::Entity& tria) {
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << tria.RefEl());
  // Fetch constant gradients of barycentric coordinate functions and area
  auto [grad_bary_coords, area] = getGradBaryCoords(tria);
  // Get values of vector field in midpoints of edges
  const Eigen::MatrixXd refqrnodes{
      (Eigen::MatrixXd(2, 3) << 0.5, 0.5, 0.0, 0.0, 0.5, 0.5).finished()};
  std::vector<Eigen::Vector2d> av_values{av_(tria, refqrnodes)};
  Eigen::Matrix<double, 1, 3> elmat{Eigen::Matrix<double, 1, 3>::Zero()};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      elmat(0, i) += (grad_bary_coords.col(i)).transpose()*av_values[j];
    }
  }
  return (area / 3) * elmat;
}

double testCDBLF(std::shared_ptr<lf::mesh::Mesh> mesh_p, Eigen::Vector2d a,
                 Eigen::Vector2d b) {
  // Build linear mesh functions
  auto vec_a = [a](Eigen::Vector2d /*x*/) -> Eigen::Vector2d { return a; };
  auto lin_b = [b](Eigen::Vector2d x) -> double { return b.dot(x); };
  const lf::mesh::utils::MeshFunctionGlobal MF_va(vec_a);
  const lf::mesh::utils::MeshFunctionGlobal MF_b(lin_b);
  // DofHandler for linear Lagrangian FE space
  lf::assemble::UniformFEDofHandler dh_linfe(mesh_p,
                                             {{lf::base::RefEl::kPoint(), 1}});
  // DofHandler for piecewise constants
  lf::assemble::UniformFEDofHandler dh_constfe(mesh_p,
                                               {{lf::base::RefEl::kTria(), 1}});
  // Matrix in triplet format holding temporary Galerkin matrix
  const unsigned int N_linfe = dh_linfe.NumDofs();
  const unsigned int N_constfe = dh_constfe.NumDofs();
  lf::assemble::COOMatrix<double> mat(N_constfe, N_linfe);
  // Initialize Galerkin matrix
  CDBLFElemMatProvider cdblfelmatprovider(MF_va);
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dh_linfe, dh_constfe, cdblfelmatprovider);

  // Compute coefficent vector of nodal interpolant
  lf::uscalfe::FeSpaceLagrangeO1<double> fes_lin(mesh_p);
  auto cv_interp = lf::fe::NodalProjection(fes_lin, MF_b);

  // Indirectly compute integral of inner product of the two vectors over domain
  const Eigen::SparseMatrix<double> A(mat.makeSparse());
  return Eigen::VectorXd::Constant(N_constfe, 1.0).transpose() * A * cv_interp;
}

}  // namespace cblfdemo

int main(int /*argc*/, char** /*argv*/) {
  // Obtain a purely triangular mesh from the collection of LehrFEM++'s
  // built-in meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};
  // vectors for testing
  const Eigen::Vector2d a(1.0, 2.0);
  const Eigen::Vector2d b(3.0, 2.0);
  // Run test computation
  double itg = cblfdemo::testCDBLF(mesh_p,a,b);
  std::cout << "Integral = " << itg << std::endl; 
  return 0;
}
