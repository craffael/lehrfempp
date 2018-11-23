#ifndef LF_FETEST_H
#define LF_FETEST_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions for testing finite element facilities
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/refinement/mesh_hierarchy.h>
#include "fe_tools.h"

namespace lf::fe {
/**
 * @brief record interpolation errors in L2 norm and H1 norm on a sequence of
 * 2D hybrid meshes
 *
 * @tparam FFUNC functor type providing the scalar-valued function to be
 *               interpolated
 * @tparam GRADFUNC functor type for objects returning the gradient
 *
 * @param mesh_ptrs array of pointers to hybrid 2D meshes
 * @param f object encoding the function
 * @param grad_f gradient of f
 * @param rfs_tria_p pointer to description of local FE space on triangles
 * @param rfs_quad_p pointer to description of local FE space on quadrilaterals
 *
 * ### type requirement
 *
 * - FFUNC must provide an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a scalar
 * - GRADFUNC must have an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a vector.
 */
template <typename FFUNC, typename GRADFUNC>
std::vector<std::pair<double, double>> InterpolationErrors(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC f,
    GRADFUNC grad_f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  // Vector of error norms
  std::vector<std::pair<double, double>> err_norms{};

  // Helper class for L2 error computation
  LocalL2NormDifference lc_L2(rfs_tria_p, rfs_quad_p, f);
  // Helper class for H1 semi norm
  LocL2GradientFEDifference lc_H1(rfs_tria_p, rfs_quad_p, grad_f);

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
    // Perform (nodal) projection of the passed function onto the finite element
    // space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, f, DefaultEntitySelector())};
    // Compute norms of interpolation error by means of numerical quadrature
    // whose order is controlled by the polynomials degree of the FE space
    double L2err = NormOfDifference(dofh, lc_L2, coeff_vec);
    double H1serr = NormOfDifference(dofh, lc_H1, coeff_vec);
    err_norms.emplace_back(L2err, H1serr);
  }
  return err_norms;
}

template <typename FFUNC, typename GRADFUNC>
inline std::vector<std::pair<double, double>> InterpolationErrors(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC f, GRADFUNC grad_f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  // Set up array of pointers to the meshes contained in the mesh hierarchy
  // std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs{};
  // // Number of levels in the hierarchy
  // lf::assemble::size_type L = multi_mesh.NumLevels();
  // for (int level = 0; level < L; level++) {
  //   // Retrieve pointer to mesh on a particular level
  //   mesh_ptrs.push_back(multi_mesh.getMesh(level));
  // }
  return InterpolationErrors(multi_mesh.getMeshes(), f, grad_f, rfs_tria_p, rfs_quad_p);
}

}  // namespace lf::fe

#endif
