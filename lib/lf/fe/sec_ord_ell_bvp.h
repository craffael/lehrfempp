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

#include "loc_comp_ellbvp.h"
#include "fe_space.h"

namespace lf::fe {
/**
 * @brief Interface class to the data defining a second-order elliptic boundary
 *        value problem in two dimensions
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
   * @brief selectors for boundary/interface conditions
   *
   * All boundary edges are visited.
   * - If EssentialConditionsOnEdge() is `true`, then the edge is treated as an
   *   an edge on the Dirichlet part of the boundary
   * - else if IsNeumannEdge() is `true`, the edge is treated as a part of the
   *   boundary where Neumann boundary conditions are imposed.
   * - else the edge is supposed to carry impedance boundary conditions
   */
  /** @{ */
  /** @brief Predicate telling whether solution is fixed on a particular edge */
  virtual bool EssentialConditionsOnEdge(
      const lf::mesh::Entity& edge) const = 0;
  /** @brief Predicate telling whether an edge carries Neumann boundary
   * conditions */
  virtual bool IsNeumannEdge(const lf::mesh::Entity& edge) const = 0;
  /** @} */

  /** @brief Data functions */
  /** @{ */
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
  bool IsNeumannEdge(const lf::mesh::Entity& /*edge*/) const override {
    return true;
  }
  double f(Eigen::Vector2d x) const override { return f_(x); }
  double g(Eigen::Vector2d /*x*/) const override { return 0.0; }
  double h(Eigen::Vector2d x) const override { return h_(x); }

 private:
  FUNCTOR_F f_;
  FUNCTOR_H h_;
};

  /** @brief Different locations for global shape functions */
enum ShapeFnType : unsigned int { kDirGSF, kNeuGSF, kImpGSF, kIntGSF };

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
 */
template <typename SCALAR>
std::pair<Eigen::SparseMatrix<SCALAR>, Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>>
SecOrdEllBVPLagrFELinSys(
    const UniformScalarFiniteElementSpace& fe_space,
    std::shared_ptr<const SecondOrderEllipticBVP<SCALAR>> bvp_p) {
  // The underlying finite element mesh
  const lf::mesh::Mesh& mesh{*fe_space.Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler& dofh{fe_space.LocGlobMap()};

  // Classify global shape functions according to their assciations with
  // different types of edges; default classification "interior"
  
  
}

}  // namespace lf::fe

#endif
