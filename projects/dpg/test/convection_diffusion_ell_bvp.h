#ifndef PROJECTS_DPG_CONVECTION_DIFFUSION_ELL_BVP
#define PROJECTS_DPG_CONVECTION_DIFFUSION_ELL_BVP

/**
 * @file
 * @brief Data structures/functions for the dpg discretization
 * of the convection diffusion problem in 2D
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../product_fe_space.h"

namespace projects::dpg::test {

/**
 * @brief Interface class to the data defining a convection-diffusion
 * boundary value problem in two dimensions.
 *
 * @tparam SCALAR a scalar type, like double
 */
template <typename SCALAR>
class ConvectionDiffusionEllipticBVP {
 protected:
  ConvectionDiffusionEllipticBVP() = default;
  ConvectionDiffusionEllipticBVP(const ConvectionDiffusionEllipticBVP&) =
      default;
  ConvectionDiffusionEllipticBVP(ConvectionDiffusionEllipticBVP&&) noexcept =
      default;
  ConvectionDiffusionEllipticBVP& operator=(
      const ConvectionDiffusionEllipticBVP&) = default;
  ConvectionDiffusionEllipticBVP& operator=(
      ConvectionDiffusionEllipticBVP&&) noexcept = default;

 public:
  /** @brief coefficient functions*/
  /** @brief diffusion coefficient (scalar valued!!) */
  [[nodiscard]] virtual SCALAR alpha(Eigen::Vector2d x) const = 0;
  /** @brief convection field coefficient */
  [[nodiscard]] virtual Eigen::Matrix<SCALAR, 2, 1> beta(
      Eigen::Vector2d x) const = 0;

  /**
   * @brief selectors for boundary conditions
   *
   * All boundary edges are visited.
   *
   * - If EssentialConditionsOnEdge() is true, Dirichlet boundary conditions are
   * imposed on that edge.
   * -  If EssentialConditionsOnEdge() is false, Neumann boundary conditions are
   * imposed.
   *
   */
  [[nodiscard]] virtual bool DirichletConditionsOnEdge(
      const lf::mesh::Entity& edge) const = 0;

  /** @brief data functions */
  /** @brief rhs source function */
  [[nodiscard]] virtual SCALAR f(Eigen::Vector2d x) const = 0;
  /** @brief dirichlet data */
  [[nodiscard]] virtual SCALAR g(Eigen::Vector2d x) const = 0;
  /** @brief neumann data */
  [[nodiscard]] virtual SCALAR h(Eigen::Vector2d x) const = 0;

  // virtual destructor:
  virtual ~ConvectionDiffusionEllipticBVP() = default;
};

/**
 * @brief class describing a full Convection diffusion boundary value problem.
 *
 * @tparam ALPHA_FUNCTOR type for passing the diffusion coefficient
 * @tparam BETA_FUNCTOR type for passing the convection field
 * @tparam F_FUNCTOR type for passing the source function
 * @tparam G_FUNCTOR type for passing the dirichlet data
 * @tparam H_FUNCTOR type for passign the neumann data
 * @tparam DIRICHLET_SELECTOR tpye to select the dirichlet edges
 */
template <typename ALPHA_FUNCTOR, typename BETA_FUNCTOR, typename F_FUNCTOR,
          typename G_FUNCTOR, typename H_FUNCTOR, typename DIRICHLET_SELECTOR>
class FullConvectionDiffusionBVP
    : public ConvectionDiffusionEllipticBVP<double> {
 public:
  /** @brief Main Constructor, constructs the problem from the
   * provided coefficients and source functions */
  FullConvectionDiffusionBVP(ALPHA_FUNCTOR alpha, BETA_FUNCTOR beta,
                             F_FUNCTOR f, G_FUNCTOR g, H_FUNCTOR h,
                             DIRICHLET_SELECTOR dirichlet_boundary)
      : alpha_(alpha),
        beta_(beta),
        f_(f),
        g_(g),
        h_(h),
        dirichlet_boundary_(std::move(dirichlet_boundary)) {}

  [[nodiscard]] double alpha(Eigen::Vector2d x) const override {
    return alpha_(x);
  }

  [[nodiscard]] Eigen::Matrix<double, 2, 1> beta(
      Eigen::Vector2d x) const override {
    return beta_(x);
  }
  [[nodiscard]] bool DirichletConditionsOnEdge(
      const lf::mesh::Entity& edge) const override {
    return dirichlet_boundary_(edge);
  }

  [[nodiscard]] double f(Eigen::Vector2d x) const override { return f_(x); }
  [[nodiscard]] double g(Eigen::Vector2d x) const override { return g_(x); }
  [[nodiscard]] double h(Eigen::Vector2d x) const override { return h_(x); }

 private:
  ALPHA_FUNCTOR alpha_;
  BETA_FUNCTOR beta_;
  F_FUNCTOR f_;
  G_FUNCTOR g_;
  H_FUNCTOR h_;
  DIRICHLET_SELECTOR dirichlet_boundary_;
};

/**
 * @brief Builds the finite element LSE for a second-order convection-diffusion
 * boundary value problem based on a DPG discretization using uniform finite
 * elements.
 *
 * @tparam SCALAR sclar type, like double
 *
 * @param fe_space the trial product fe-space of the  DPG method.
 * @param element_matrix_provider provides the DPG element matrix
 * @param element_vector_provider provides the DPG element vector
 * @param bvp_p pointer to a description of the convection-diffusion boundary
 * value problem
 * @param trace_component index of the trace-like component on which the
 * Dirichlet boundary conditions are imposed (as essential boundary conditions)
 * @param flux_component index of the flux-like component on which the Neumann
 * boundary conditions are imposed (as essential boundary conditions)
 *
 * @return A pair containing the sparse Galerkin matix and the right hand side
 * vector
 *
 */
template <typename SCALAR>
std::pair<Eigen::SparseMatrix<SCALAR>, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>
ConvectionDiffusionDPGLinSys(
    std::shared_ptr<ProductUniformFESpace<SCALAR>> fe_space,
    std::shared_ptr<DpgElementMatrixProvider<SCALAR>> element_matrix_provider,
    std::shared_ptr<DpgElementVectorProvider<SCALAR>> element_vector_provider,
    std::shared_ptr<const ConvectionDiffusionEllipticBVP<SCALAR>> bvp_p,
    size_type trace_component, size_type flux_component) {
  // extract mesh and dofhandler:
  const lf::mesh::Mesh& mesh{*fe_space->Mesh()};
  const ProductUniformFEDofHandler& dofh{fe_space->LocGlobMap()};

  // assemble the full galerkin matrix
  const size_type N_dofs(dofh.NumDofs());
  lf::assemble::COOMatrix<SCALAR> A(N_dofs, N_dofs);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, *element_matrix_provider,
                                      A);

  // std::cout << "Assembled A " << std::endl;
  // assemble right hand side:
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  lf::assemble::AssembleVectorLocally(0, dofh, *element_vector_provider, phi);
  // std::cout << "Assembled phi" << std::endl;

  // assemble boundary conditions:
  // obtain pointers to reference shape functions on edges:
  auto rfs_edge_trace_p = fe_space->ShapeFunctionLayout(
      lf::base::RefEl::kSegment(), trace_component);
  auto rfs_edge_flux_p = fe_space->ShapeFunctionLayout(
      lf::base::RefEl::kSegment(), flux_component);
  LF_ASSERT_MSG(rfs_edge_trace_p != nullptr && rfs_edge_flux_p != nullptr,
                "missing rfs on edges");

  // obtain selectors for dirichlet and neumann data on the boundary:
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  auto dirichlet_selector = [&bd_flags,
                             &bvp_p](const lf::mesh::Entity& edge) -> bool {
    return bd_flags(edge) && bvp_p->DirichletConditionsOnEdge(edge);
  };
  auto neumann_selector = [&bd_flags,
                           &bvp_p](const lf::mesh::Entity& edge) -> bool {
    return bd_flags(edge) && !bvp_p->DirichletConditionsOnEdge(edge);
  };

  // count the number of dirichlet and neumann edges:
  size_type no_Dirichlet_edges = 0;
  size_type no_Neumann_edges = 0;
  for (const lf::mesh::Entity* const edge : mesh.Entities(1)) {
    if (dirichlet_selector(*edge)) {
      no_Dirichlet_edges++;
    }
    if (neumann_selector(*edge)) {
      no_Neumann_edges++;
    }
  }

  // std::cout << "Assembly with Dirichlet edges: " << no_Dirichlet_edges << ",
  // Neumann edges: " <<
  //           no_Neumann_edges << std::endl;

  // wrap the boundary conditions into mesh functions
  auto g_mf = lf::mesh::utils::MeshFunctionGlobal(
      [&bvp_p](auto x) -> SCALAR { return bvp_p->g(x); });
  auto h_mf = lf::mesh::utils::MeshFunctionGlobal(
      [&bvp_p](auto x) -> SCALAR { return bvp_p->h(x); });

  // obtain flags and values of degrees of feedom.
  auto ess_bdc_flags_values{InitEssentialConditionsFromFunctions(
      dofh, *rfs_edge_trace_p, *rfs_edge_flux_p, dirichlet_selector,
      neumann_selector, g_mf, h_mf, trace_component, flux_component)};

  // eliminate dirichlet and neumann dofs from the linear system:
  lf::assemble::FixFlaggedSolutionComponents<SCALAR>(
      [&ess_bdc_flags_values](lf::base::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A, phi);

  return {A.makeSparse(), phi};
}

}  // namespace projects::dpg::test

#endif  // PROJECTS_DPG_CONVECTION_DIFFUSION_ELL_BVP
