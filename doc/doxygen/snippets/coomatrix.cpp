/**
 * @file
 * @brief Doxygen snippets to show lf::assemble::COOMatrix usage
 * @author Ralf Hiptmair
 * @date Tue 29 Oct 2019 03:09:10 PM CET
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::assemble {
void snippetfoo() {
  //! [usage]
  // Obtain a shared pointed to a mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0)};
  // Initialization of index mapping for linear finite elements
  lf::assemble::UniformFEDofHandler loc_glob_map(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});
  // Initialize ENTITYMATRIXPROVIDER object for local computations
  lf::uscalfe::LinearFELaplaceElementMatrix loc_mat_laplace{};
  // Dimension of finite element space
  const lf::assemble::size_type N_dofs(loc_glob_map.NoDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);
  // Building the Galerkin matrix (trial space = test space)
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, loc_glob_map, loc_glob_map, loc_mat_laplace, mat);
  // the `mat` object now contains the Galerkin matrix in triplet/COO format
  // Initialize an Eigen sparse matrix object from `mat`
  Eigen::SparseMatrix<double> stiffness_matrix = {mat.makeSparse()};
  //! [usage]
}

}  // namespace lf::assemble
