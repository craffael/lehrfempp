/**
 * @file
 * @brief Computation of L2-norm of piecewise linear FE function
 * assemble module; meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include "lecturedemotwonorm.h"

namespace lecturedemo {

class LinFEMassMatrixProvider {
 public:
  /** @brief default constructor */
  explicit LinFEMassMatrixProvider() = default;
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /** @brief Main method for computing the element vector
   * @param cell refers to current cell for which the element vector is desired
   * The implementation uses an analytic formula defined over triangular cells*/
  Eigen::Matrix3d Eval(const lf::mesh::Entity &tria);

 private:
  // clang-format off
  const Eigen::Matrix3d elMat_{(Eigen::Matrix3d() <<
				   2.0, 1.0, 1.0,
				   1.0, 2.0, 1.0,
				   1.0, 1.0, 2.0).finished()};

  // clang-format on
};

/** Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix3d LinFEMassMatrixProvider::Eval(const lf::mesh::Entity &tria) {
  // Throw error in case no triangular cell
  LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << tria.RefEl());
  // Compute the area of the triangle cell
  const double area = lf::geometry::Volume(*(tria.Geometry()));
  // Compute the mass element matrix over the cell
  return (area / 12.0 * elMat_);  // return the local mass element matrix
}

double l2normByMassMatrix(const lf::assemble::DofHandler &dofh,
                          const Eigen::VectorXd &uvec) {
  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NoDofs());
  LF_ASSERT_MSG(
      N_dofs == uvec.size(),
      "Size mismatch: NoDofs = " << N_dofs << " <-> size = " << uvec.size());
  // Obtain Galerkin mass matrix by local cell-oriented assembly
  LinFEMassMatrixProvider M_loc{};  // ENTITY_MATRIX_PROVIDER
  lf::assemble::COOMatrix<double> M_coo{
      lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
          0, dofh, M_loc)};
  // Optional: output of mass matrix for debugging purposes
  // std::cout << "Mass matrix" << std::endl << M_coo.makeDense() << std::endl;
  // Multiply coefficient vector onto mass matrix from left and right
  return std::sqrt(uvec.dot(M_coo.MatVecMult(1.0, uvec)));
}

double l2normByQuadrature(const lf::assemble::DofHandler &dofh,
                          const Eigen::VectorXd &uvec) {
  LF_ASSERT_MSG(dofh.NoDofs() == uvec.size(),
                "Size mismatch: NoDofs = " << dofh.NoDofs()
                                           << " <-> size = " << uvec.size());
  // Obtain reference to current mesh
  const lf::mesh::Mesh &mesh{*dofh.Mesh()};
  // Summation variable for square of L2-norm
  double l2n_square{0.0};
  // Loop over the cells
  for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
    // Obtain cell type
    const lf::base::RefEl ref_el{cell.RefEl()};
    // Only lowest-order Lagrangian FE are supported
    LF_ASSERT_MSG(dofh.NoLocalDofs(cell) == ref_el.NumNodes(),
                  "No nodes must be equal to number of local shape functions");
    // Obtain shape of cell
    const lf::geometry::Geometry *geo_p{cell.Geometry()};
    // Query volume of the cell
    const double area = lf::geometry::Volume(*geo_p);
    // Obtain global indices of global shape functions covering the cell
    nonstd::span<const lf::assemble::gdof_idx_t> idx{
        dofh.GlobalDofIndices(cell)};
    switch (ref_el) {
      case lf::base::RefEl::kTria(): {
        // Edge-midpoint based local quadrature, exact for quadratic
        // polynomials, hence exact for the sqaure of piecewise linear finite
        // element functions
        const Eigen::Vector3d uloc(uvec[idx[0]], uvec[idx[1]], uvec[idx[2]]);
        const Eigen::Vector3d mpv(uloc[0] + uloc[1], uloc[1] + uloc[2],
                                  uloc[2] + uloc[0]);
        // Factor 1/12 = 1/3*1/4, 1/3|K| from quadrature weight, 1/4 from
        // squaring the factor 1/2 occurring in the midpoint value
        l2n_square += 1 / 12.0 * area * mpv.squaredNorm();
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Unsupport cell type " << ref_el);
        break;
      }
    }  // end switch
  }
  return std::sqrt(l2n_square);
}

double l2normByMeshFunction(
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
    const Eigen::VectorXd &uvec) {
  // Compute norms of finite element function by means of numerical quadrature
  // which is exact for piecewise quadratic polynomials (degree of exactness 2)
  auto mf_fe = lf::uscalfe::MeshFunctionFE<double, double>(fe_space, uvec);
  return std::sqrt(lf::uscalfe::IntegrateMeshFunction(*fe_space->Mesh(),
                                                      squaredNorm(mf_fe), 2));
}

void lecturedemotwonorm() {
  // Debugging flags for AssembleMatrixLocally()
  // lf::assemble::ass_mat_dbg_ctrl = 31;
  std::cout << "LehrFEM++ demo: computation of L2-norm of piecewise linear FE "
               "function"
            << std::endl;
  // Obtain a purely triangular mesh of the unit square from the collection of
  // LehrFEM++'s built-in meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0)};
  // Optional: Output information on the mesh
  // std::cout << "Mesh: " << std::endl << *mesh_p << std::endl;

  // Build a lowest-order Lagrangian finite element space
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Fetch the d.o.f. handler from the finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  // Optional: Output information on the d.o.f. handler
  // std::cout << "DofHandler" << std::endl << dofh << std::endl;

  // Dimension of  the finite element space
  size_type n_dofs{dofh.NoDofs()};

  // Build finite element coefficient vector by interpolating
  // a known linear (!) function. The L2-norm of linear functions
  // should be computed exactly.
  auto u = lf::uscalfe::MeshFunctionGlobal(
      [](auto x) -> double { return 2 * x[0] + x[1]; });
  const Eigen::VectorXd uvec =
      lf::uscalfe::NodalProjection<double>(*fe_space_p, u);
  // Other, simpler, ways to set the coefficient vector
  // const Eigen::VectorXd uvec{Eigen::VectorXd::LinSpaced(n_dofs,0.0,1.0)};
  // const Eigen::VectorXd uvec{Eigen::VectorXd::Random(n_dofs,1.0)};

  // Console output of different functions: results should all be equal
  // up to machine precision.
  std::cout << "Pw. linear FE function from " << n_dofs
            << "-dimenmsional FE space" << std::endl;
  std::cout << "Euclidean norm = " << uvec.norm() << std::endl;

  std::cout << "l2normByMassMatrix = " << l2normByMassMatrix(dofh, uvec)
            << std::endl;
  std::cout << "l2normByQuadrature = " << l2normByQuadrature(dofh, uvec)
            << std::endl;
  std::cout << "l2normByMeshFunction = "
            << l2normByMeshFunction(fe_space_p, uvec) << std::endl;
}

}  // namespace lecturedemo
