#ifndef HLDO_SPHERE_DEBUG_DIRAC_CONVERGENCE_H
#define HLDO_SPHERE_DEBUG_DIRAC_CONVERGENCE_H

/**
 * @file dirac_convergence_test.h
 *
 * @brief test the convergence of the solutions for the dirac operator of they
 * convergence but not to what funciton.
 */
#include <dirac_operator_source_problem.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <mesh_function_whitney_one.h>
#include <mesh_function_whitney_two.h>
#include <mesh_function_whitney_zero.h>
#include <norms.h>
#include <results_processing.h>
#include <sphere_triag_mesh_builder.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace projects::hldo_sphere {
namespace debugging {

using complex = std::complex<double>;

/**
 * @brief stores solutions for a number of refinement levels
 */
struct SolutionList {
  std::vector<Eigen::Matrix<complex, Eigen::Dynamic, 1>> mu_zero;
  std::vector<Eigen::Matrix<complex, Eigen::Dynamic, 1>> mu_one;
  std::vector<Eigen::Matrix<complex, Eigen::Dynamic, 1>> mu_two;
};

/**
 * @brief Class to test the convergence of the Dirac operator
 *
 * This test does only test convergence but not the convergence to a specific
 * value. Morespecifically its used to test
 * @f[
 *  |\|u_k\| - \|u_{k-1}\|| \to 0
 *  @f]
 *
 * The method Compute() computes a list of values stored in
 * ./results/result_dirac_convergence_`refimement_level`.csv that can be used
 * for plots
 *
 */
class DiracConvergenceTest {
 public:
  /**
   *
   * @brief Constructor setting all the functions and the reference k
   *
   * @param f_zero load function corresponding to the analytical solution u_zero
   * @param f_one load function corresponding to the analytical solution u_one
   * @param f_two load function corresponding to the analytical solution u_two
   * @param k reference used in all the functions such that changes of k affect
   * the functions
   *
   */
  DiracConvergenceTest(
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_zero,
      std::function<
          Eigen::Matrix<complex, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f_one,
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_two,
      double &k)
      : f_zero_(f_zero), f_one_(f_one), f_two_(f_two), k_(k) {}

  /**
   *
   * @brief Solves the dirac opeartor source problems up to the given
   * refinement_level
   *
   * @param refinement_levels integer maximal refinement level to compute
   *
   */
  void Compute(unsigned refinement_levels);

 private:
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_zero_;
  std::function<Eigen::Matrix<complex, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f_one_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_two_;
  double &k_;
};
}  // namespace debugging
}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_DEBUG_DIRAC_CONVERGENCE_H
