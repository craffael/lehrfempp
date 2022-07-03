#ifndef HLDO_SPHERE_OPERATORS_DIRAC_OPERATOR_SOURCE_PROBLEM_H
#define HLDO_SPHERE_OPERATORS_DIRAC_OPERATOR_SOURCE_PROBLEM_H

/**
 * @file dirac_operator_source_problem.h
 *
 * @brief Class to discretize the Dirac operator
 */

#include <dirac_operator.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh_interface.h>
#include <lf/quad/quad.h>
#include <load_vector_provider.h>
#include <mass_matrix_provider.h>
#include <sphere_triag_mesh_builder.h>
#include <whitney_one_curl_curl_matrix_provider.h>
#include <whitney_one_grad_matrix_provider.h>
#include <whitney_one_mass_matrix_provider.h>
#include <whitney_two_mass_matrix_provider.h>

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace operators {

using complex = std::complex<double>;

/**
 * @brief Computes the Galerkin LSE for the Dirac Operator source problem
 *
 * @f[
 *  D \vec{u} + \imath k \vec{u} = \vec{f}
 * @f]
 *
 * The computation must be complex valued since the problem definition involves
 * a complex @f$\imath@f$
 *
 * @note Only triangular meshes are supported
 *
 */
class DiracOperatorSourceProblem {
 public:
  /**
   * @brief Constructor
   * creates basic mesh (Octaeder with radius 1.0)
   *
   * set k to 1
   *
   * For all three Hodge Laplacians
   * creates zerovalued functions f
   *
   */
  DiracOperatorSourceProblem();

  /**
   * @brief Computes the Galerkin LSE for the dirac operator source problem
   *
   * @f[
   *   D \vec{u} + \imath k \vec{u} = \vec{f}
   * @f]
   *
   * The Galerkin matrix will be accessable with GetGalerkinMatrix()
   * The load vector will be accessable with GetLoadVector()
   *
   */
  void Compute();

  /**
   *
   * @brief solves the linear systems build in Compute
   *
   * Uses the Matrices stored in coo_matrix_ and
   * the righthandside vectors stored in phi_
   * to compute the basis expansion
   * coefficiants mu_ approximating
   *
   * @note the method Compute() has to be called prior to
   * calling Solve()
   *
   * @f[
   *   D \vec{u} + \imath k \vec{u} = \vec{f}
   * @f]
   *
   */
  void Solve();

  /**
   * @brief Sets the mesh
   * @param mesh_p pointer to the mesh
   *
   * requries all cells in the mesh are triangles
   * requries mesh global dimension to be 3
   */
  void SetMesh(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

  /**
   * @brief Sets the load functions
   * @param f0 load function in @f$ L^2 @f$
   * @param f1 load functions in @f$ L^2_t @f$ tangential vector field on the
   * sphere
   * @param f2 load functions in @f$ L^2 @f$
   */
  void SetLoadFunctions(
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f0,
      std::function<
          Eigen::Matrix<complex, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f1,
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f2) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
  }

  /**
   * @brief Sets k for the souce problems
   *
   * @param k real valued parameter scaling the mass matrix
   */
  void SetK(double k) { k_ = std::complex<double>(k, 0); }

  /**
   * @brief returns the Loadvector
   *
   * This is the righthandside of LSE
   *
   * @note The loadvector must be computed with `Compute` before calling
   * this function
   *
   */
  Eigen::Matrix<complex, Eigen::Dynamic, 1> GetLoadVector() { return phi_; }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE
   *
   * @note The Galerkin matrix must be computed with `Compute` before
   * calling this funciton
   *
   */
  lf::assemble::COOMatrix<complex> GetGalerkinMatrix() { return coo_matrix_; }

  /**
   * @brief retunrs the basis expansion coefficiants for
   * the solution of the source problem
   *
   * BUT ONLY THE REAL PART
   *
   * The basis expansion coefficiants for the first
   * number of vertices entries are with respect to the
   * barycentric basis functions on the mesh.
   *
   * The basis expansion coefficiants for the next
   * number of edge entries are with respect to the
   * whitney one form basis functions.
   *
   * The last number of cells basis expansion coefficiants
   * are with respect to the cellwise constant basis functions
   *
   * @param index of the basis (0 -> barycentric basis, 1 -> Whitney 1-forms,
   * surface edge elements, 2 -> cellwise constant)
   *
   * @returns basis expansion coefficiants
   *
   * @note the methods Compute() and Solve() must be
   * called before the result is available
   *
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> GetMu(int index);

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f0_;
  std::function<Eigen::Matrix<complex, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f1_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f2_;
  lf::assemble::COOMatrix<complex> coo_matrix_;
  Eigen::Matrix<complex, Eigen::Dynamic, 1> phi_;
  Eigen::Matrix<complex, Eigen::Dynamic, 1> mu_;
  complex k_;
};

}  // namespace operators

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_OPERATORS_DIRAC_OPERATOR_SOURCE_PROBLEM_H
