#ifndef LF_USCALFE_PROLONGATION_H
#define LF_USCALFE_PROLONGATION_H

#include <lf/assemble/dofhandler.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/mesh_function_transfer.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/mesh/utils/utils.h>

namespace lf::uscalfe {

/**
 * @brief Interpolate a vector of DOFs on a finer mesh
 * @tparam SCALAR_COEFF The scalar of the coefficient vector
 * @tparam FES_COARSE The FE space on the coarse mesh
 * @tparam FES_FINE The FE space on the fine mesh
 * @param mh A reference to the MeshHierarchy containing the underlying meshes
 * @param fespace_coarse The FE space on the coarse mesh
 * @param fespace_fine The FE space on the fine mesh
 * @param dofs_coarse A basis function coefficient vector on the coarse mesh
 * @param level The level of the coarse mesh
 * @returns An interpolated vector of basis function coefficients on the fine
 * mesh
 */
template <typename SCALAR_COEFF, typename FES_COARSE, typename FES_FINE>
[[nodiscard]] Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1> prolongate(
    const lf::refinement::MeshHierarchy &mh,
    std::shared_ptr<const FES_COARSE> fespace_coarse,
    std::shared_ptr<const FES_FINE> fespace_fine,
    const Eigen::Matrix<SCALAR_COEFF, Eigen::Dynamic, 1> &dofs_coarse,
    lf::base::size_type level_coarse) {
  // Assert that the FES_* are actually FE spaces
  using scalar_fe_coarse_t = typename FES_COARSE::Scalar;
  using scalar_fe_fine_t = typename FES_FINE::Scalar;
  static_assert(
      std::is_convertible_v<
          FES_COARSE, lf::uscalfe::UniformScalarFESpace<scalar_fe_coarse_t>>,
      "Invalid coarse FE space provided");
  static_assert(
      std::is_convertible_v<
          FES_FINE, lf::uscalfe::UniformScalarFESpace<scalar_fe_fine_t>>,
      "Invalid fine FE space provided");
  // Obtain the dofhandlers from the fe spaces
  const lf::assemble::DofHandler &dofh_coarse{fespace_coarse->LocGlobMap()};
  const lf::assemble::DofHandler &dofh_fine{fespace_fine->LocGlobMap()};
  const lf::base::size_type N_coarse = dofh_coarse.NumDofs();
  const lf::base::size_type N_fine = dofh_fine.NumDofs();
  // Assert correctness of inputs
  LF_ASSERT_MSG(level < mh.NumLevels() - 1,
                "level must not point to the finest mesh in the hierarchy");
  LF_ASSERT_MSG(
      dofs_coarse.size() >= N_coarse,
      "Too few basis function coefficients provided for coarse FE space");
  // Construct a mesh function to simplify the point evaluations
  const lf::uscalfe::MeshFunctionFE mf_coarse(fespace_coarse, dofs_coarse);
  // Transfer the mesh function to the finer mesh
  const lf::refinement::MeshFunctionTransfer mf_fine(mh, mf_coarse, level_coarse);
  // Return the nodal projection of this transferred mesh function
  return lf::uscalfe::NodalProjection(*fespace_fine, mf_fine);
}

}   // end namespace lf::uscalfe

#endif // LF_USCALFE_PROLONGATION_H
