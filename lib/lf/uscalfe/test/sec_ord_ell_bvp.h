#ifndef LF_ELLBVP_H
#define LF_ELLBVP_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structures/functions for the Lagrangian finite element
 * discretization of second order elliptic boundary value problems in 2D
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/lagrfe.h>

namespace lf::uscalfe::test {
/**
 * @brief Interface class to the data defining a second-order elliptic boundary
 *        value problem in two dimensions
 *
 * @tparam SCALAR a scalar type, usually `double` or `complex<double>`
 */
template <typename SCALAR>
class SecondOrderEllipticBVP {
 protected:
  SecondOrderEllipticBVP() = default;
  SecondOrderEllipticBVP(const SecondOrderEllipticBVP&) = default;
  SecondOrderEllipticBVP(SecondOrderEllipticBVP&&) noexcept = default;
  SecondOrderEllipticBVP& operator=(const SecondOrderEllipticBVP&) = default;
  SecondOrderEllipticBVP& operator=(SecondOrderEllipticBVP&&) noexcept =
      default;

 public:
  /** @brief coefficient functions */
  /** @{ */
  /** @brief diffusion coefficient tensor */
  virtual Eigen::Matrix<SCALAR, 2, 2> alpha(Eigen::Vector2d x) const = 0;
  /** @brief reaction coefficient  */
  virtual SCALAR gamma(Eigen::Vector2d x) const = 0;
  /** @brief impedance coefficient */
  virtual SCALAR eta(Eigen::Vector2d x) const = 0;
  /** @} */

  /**
   * @defgroup sel
   * @brief selectors for boundary/interface conditions
   *
   * All boundary edges are visited.
   * - If EssentialConditionsOnEdge() is `true`, then the edge is treated as an
   *   an edge on the Dirichlet part of the boundary
   * - else if IsImpedanceEdge() is `true`, the edge is treated as a part of
   * the boundary where impedance boundary conditions are imposed.
   * - else the edge is supposed to belong to the Neumann boundary part
   * @{ */
  /** @brief Predicate telling whether solution is fixed on a particular edge */
  virtual bool EssentialConditionsOnEdge(
      const lf::mesh::Entity& edge) const = 0;
  /** @brief Predicate telling whether an edge carries Impendance boundary
   * conditions */
  virtual bool IsImpedanceEdge(const lf::mesh::Entity& edge) const = 0;
  /** @} */

  /** @defgroup data
   * @brief Data functions
   * @{ */
  /** @brief right-hand side source function in @f$\in L^2 @f$ */
  virtual SCALAR f(Eigen::Vector2d x) const = 0;
  /** @brief Dirichlet data @f$\in C^0 @f$ */
  virtual SCALAR g(Eigen::Vector2d x) const = 0;
  /** @brief Neumann and impedance data @f$\in L^2 @f$ */
  virtual SCALAR h(Eigen::Vector2d x) const = 0;
  /** @} */

  virtual ~SecondOrderEllipticBVP() = default;
};

/**
 * @brief Class describing a pure Neumann problem for the Laplacian
 *
 * @tparam FUNCTOR_F type for passing the source function
 * @tparam FUNCTOR_H type for passing the Neumann data
 *
 * Specifies a second-order scalar elliptic boundary value problem with
 * pure inhomogeneous Neumann boundary conditions.
 *
 */
template <typename FUNCTOR_F, typename FUNCTOR_H>
class PureNeumannProblemLaplacian : public SecondOrderEllipticBVP<double> {
 public:
  PureNeumannProblemLaplacian(FUNCTOR_F f, FUNCTOR_H h) : f_(f), h_(h) {}
  PureNeumannProblemLaplacian(const PureNeumannProblemLaplacian&) = default;
  PureNeumannProblemLaplacian(PureNeumannProblemLaplacian&&) noexcept = default;
  PureNeumannProblemLaplacian& operator=(const PureNeumannProblemLaplacian&) =
      default;
  PureNeumannProblemLaplacian& operator=(
      PureNeumannProblemLaplacian&&) noexcept = default;

  Eigen::Matrix<double, 2, 2> alpha(Eigen::Vector2d x) const override {
    return Eigen::Matrix<double, 2, 2>::Identity(2, 2);
  }
  double gamma(Eigen::Vector2d x) const override { return 0.0; }
  double eta(Eigen::Vector2d x) const override { return 0.0; }
  bool EssentialConditionsOnEdge(
      const lf::mesh::Entity& /*edge*/) const override {
    return false;
  }
  bool IsImpedanceEdge(const lf::mesh::Entity& /*edge*/) const override {
    return false;
  }
  double f(Eigen::Vector2d x) const override { return f_(x); }
  double g(Eigen::Vector2d /*x*/) const override { return 0.0; }
  double h(Eigen::Vector2d x) const override { return h_(x); }

 private:
  FUNCTOR_F f_;
  FUNCTOR_H h_;
};

// TODO(ralfh) Putting this into a header file leads to multiple definitions of
// the same symbol!!!
// /** @brief output control variable for function SecOrdEllBVPLagrFELinSys() */
// CONTROLDECLAREINFO(
//     LFELinSys_ctrl, "LFELinSys_ctrl",
//     "Output control variable for function SecOrdEllBVPLagrFELinSys()");
// static const unsigned int kLFELinSys_bdinfo = 2;

/**
 * @brief Builds finite element linear system of equations for a second-order
 * elliptic boundary value problem using Lagrangian finite elements of uniform
 * degree.
 *
 * @tparam SCALAR scalar type for entries of Galerkin matrix and right-hand side
 * vector
 *
 * @param fe_space Lagrangian finite element space of uniform polynomial degree
 * @param bvp_p shared pointer to description of boundary value problem
 * @return both the sparse finite element Galerkin matrix and the right-hand
 * side vector
 *
 * It relies on the various assembly functions provided in fe_tools.h
 * - lf::uscalfe::SecOrdBVPLagrFEFullInteriorGalMat()
 * - lf::uscalfe::SecOrdBVPLagrFEBoundaryGalMat()
 * - lf::uscalfe::LagrFEVolumeRightHandSideVector()
 * - lf::uscalfe::LagrFEBoundaryRightHandSideVector()
 */
template <typename SCALAR>
std::pair<Eigen::SparseMatrix<SCALAR>, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>
SecOrdEllBVPLagrFELinSys(
    std::shared_ptr<ScalarUniformFESpace<SCALAR>> fe_space,
    std::shared_ptr<const SecondOrderEllipticBVP<SCALAR>> bvp_p) {
  LF_ASSERT_MSG(bvp_p != nullptr, "No valid BVP specified");

  // The underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // Preprocessing: count number of edges with different boundary conditions
  size_type no_Dirichlet_edges = 0;
  size_type no_Neumann_edges = 0;
  size_type no_impedance_edges = 0;
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  for (const lf::mesh::Entity& edge : mesh.Entities(1)) {
    if (bd_flags(edge)) {
      if (bvp_p->EssentialConditionsOnEdge(edge)) {
        no_Dirichlet_edges++;
      } else if (bvp_p->IsImpedanceEdge(edge)) {
        no_impedance_edges++;
      } else {
        no_Neumann_edges++;
      }
    }
  }

  // TODO(ralfh): Fix this once LFELinSys_ctrl has been fixed.
  // SWITCHEDSTATEMENT(LFELinSys_ctrl, kLFELinSys_bdinfo,
  //                   std::cout << "B.c.: " << no_Dirichlet_edges
  //                             << " Dirichlet edges, " << no_Neumann_edges
  //                             << " Neumann edges, " << no_impedance_edges
  //                             << " impedance edges" << std::endl);

  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NoDofs());
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<SCALAR> A(N_dofs, N_dofs);

  // I: Assemble finite element Galerkin matrix
  // First the volume part for the bilinear form
  SecOrdBVPLagrFEFullInteriorGalMat(
      fe_space,
      MeshFunctionGlobal([&bvp_p](auto x) -> Eigen::Matrix<SCALAR, 2, 2> {
        return (bvp_p->alpha(x));
      }),
      MeshFunctionGlobal(
          [&bvp_p](auto x) -> double { return (bvp_p->gamma(x)); }),
      A);
  // Update with potential contributions from edges (Impedance boundary
  // conditions)
  if (no_impedance_edges > 0) {
    SecOrdBVPLagrFEBoundaryGalMat(
        fe_space, MeshFunctionGlobal([&bvp_p](auto x) -> SCALAR {
          return (bvp_p->eta(x));
        }),
        [&bvp_p](const lf::mesh::Entity& edge) -> bool {
          return (!bvp_p->EssentialConditionsOnEdge(edge)) &&
                 (bvp_p->IsImpedanceEdge(edge));
        },
        A);
  }

  // II: Right-hand side vector; has to be set to zero initially
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // Assemble volume part of right-hand side vector depending on the function f
  LagrFEVolumeRightHandSideVector(
      fe_space,
      MeshFunctionGlobal([&bvp_p](auto x) -> SCALAR { return (bvp_p->f(x)); }),
      phi);
  // Add contributions from Neumann and impedance edges
  if ((no_Neumann_edges > 0) || (no_impedance_edges > 0)) {
    LagrFEBoundaryRightHandSideVector(
        fe_space, MeshFunctionGlobal([&bvp_p](auto x) -> SCALAR {
          return (bvp_p->h(x));
        }),
        [&bvp_p, &bd_flags](const lf::mesh::Entity& edge) -> bool {
          return (bd_flags(edge) && (!bvp_p->EssentialConditionsOnEdge(edge)));
        },
        phi);
  }

  // III: Fixing coefficients due to essential boundary conditions
  if (no_Dirichlet_edges > 0) {
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_edge_p =
        fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
    LF_ASSERT_MSG(rfs_edge_p != nullptr, "FE specification for edges missing");

    // Obtain flags and values for degrees of freedom located on Dirichlet edges
    auto ess_bdc_flags_values{InitEssentialConditionFromFunction(
        dofh, *rfs_edge_p,
        [&bvp_p, &bd_flags](const lf::mesh::Entity& edge) -> bool {
          return (bd_flags(edge) && bvp_p->EssentialConditionsOnEdge(edge));
        },
        MeshFunctionGlobal(
            [&bvp_p](auto x) -> SCALAR { return bvp_p->g(x); }))};
    // Eliminate Dirichlet dofs from linear system
    lf::assemble::fix_flagged_solution_components<SCALAR>(
        [&ess_bdc_flags_values](glb_idx_t gdof_idx) {
          return ess_bdc_flags_values[gdof_idx];
        },
        A, phi);
  }

  return {A.makeSparse(), phi};
}

}  // namespace lf::uscalfe::test

#endif
