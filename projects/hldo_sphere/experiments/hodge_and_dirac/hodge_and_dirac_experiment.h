#ifndef THESIS_EXPERIMENTS_HODGE_AND_DIRAC_H
#define THESIS_EXPERIMENTS_HODGE_AND_DIRAC_H

/**
 * @file hodge_and_dirac_experiment.h
 */
#include <dirac_operator_source_problem.h>
#include <hodge_laplacians_source_problems.h>
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

namespace projects::hldo_sphere::experiments {

using complex = std::complex<double>;

/**
 * @brief Creates and solves the discretised Hodge Laplacian source problems
 * and the Dirac operator source Problem for
 * a given list of levels and values of k and compares the two solutions.
 *
 * In order to make this work, the load functions for the Dirac operator
 * and Hodge Laplacians must be appropriately chosen
 *
 * @f[
 *  D\, \vec{u} + \imath\, k \, \vec{u} = \vec{f} \\
 *  \Delta\, \vec{u} + k^2\, \vec{u} = D\, \vec{f} - \imath k \, \vec{f}
 * @f]
 *
 * Details about the experiments can be found in the thesis
 * `Hodge-Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * chapter 6.
 *
 */
class HodgeAndDiracExperiment {
 public:
  /**
   *
   * @brief Constructor setting all the functions and the reference k
   *
   * @param f_zero load function to the analytical solution of (D + ik) u_zero
   * @param f_one load function to the analytical solution (D + ik) u_one
   * @param f_two load function to the analytical solution (D + ik) u_two
   * @param f_zero_dirac load function equivalent to  (D - ik) f_zero
   * @param f_one_dirac load function  equivalent to (D - ik) f_one
   * @param f_two_dirac load function equivalent to (D - ik) f_two
   * @param k reference used in all the functions such that changes of k affect
   * the functions
   * @param name identifier of the example (cretes a folder with this name for
   * the results)
   *
   */
  HodgeAndDiracExperiment(
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_zero_dirac,
      std::function<
          Eigen::Matrix<complex, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f_one_dirac,
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_two_dirac,
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_zero,
      std::function<
          Eigen::Matrix<complex, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f_one,
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_two,
      double &k, std::string name)
      : f_zero_dirac_(f_zero_dirac),
        f_one_dirac_(f_one_dirac),
        f_two_dirac_(f_two_dirac),
        f_zero_(f_zero),
        f_one_(f_one),
        f_two_(f_two),
        k_(k),
        name_(name) {}

  /**
   *
   * @brief Solves the Hodge Laplacian and the Dirac operator source problems
   * for the tensor product of passed refinement levels and ks
   *
   * @param refinement_levels integer list containing all the levels
   * @param ks list of all ks to be used
   *
   */
  void Compute(std::vector<unsigned> refinement_levels,
               std::vector<double> ks) {
    int nl = refinement_levels.size();
    int nk = ks.size();

    // Initialize solution wrapper
    projects::hldo_sphere::post_processing::ProblemSolutionWrapper<complex>
        solutions_lap;
    solutions_lap.k = ks;
    solutions_lap.levels = refinement_levels;
    solutions_lap.mesh = std::vector<std::shared_ptr<const lf::mesh::Mesh>>(nl);
    solutions_lap.solutions = std::vector<
        projects::hldo_sphere::post_processing::ProblemSolution<complex>>(nl);

    projects::hldo_sphere::post_processing::ProblemSolutionWrapper<complex>
        solutions_dirac;
    solutions_dirac.k = ks;
    solutions_dirac.levels = refinement_levels;
    solutions_dirac.mesh =
        std::vector<std::shared_ptr<const lf::mesh::Mesh>>(nl);
    solutions_dirac.solutions = std::vector<
        projects::hldo_sphere::post_processing::ProblemSolution<complex>>(nl);

    projects::hldo_sphere::operators::HodgeLaplaciansSourceProblems<complex>
        hodge_builder;

    projects::hldo_sphere::operators::DiracOperatorSourceProblem dirac_builder;

    // functions are passed by reference hence changing the k still influences
    // the functions
    hodge_builder.SetLoadFunctions(f_zero_dirac_, f_one_dirac_, f_two_dirac_);
    dirac_builder.SetLoadFunctions(f_zero_, f_one_, f_two_);

    // Read the mesh from the gmsh file
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

    projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
        projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

    // we always take the same radius
    sphere.setRadius(1.);

    // start timer for total time
    std::chrono::steady_clock::time_point start_time_total =
        std::chrono::steady_clock::now();

    for (unsigned il = 0; il < nl; ++il) {
      unsigned lvl = refinement_levels[il];
      // Initialize solution vectors
      solutions_lap.solutions[il].mu_zero.resize(nk);
      solutions_lap.solutions[il].mu_one.resize(nk);
      solutions_lap.solutions[il].mu_two.resize(nk);

      solutions_dirac.solutions[il].mu_zero.resize(nk);
      solutions_dirac.solutions[il].mu_one.resize(nk);
      solutions_dirac.solutions[il].mu_two.resize(nk);

      for (int ik = 0; ik < nk; ik++) {
        k_ = ks[ik];
        std::cout << "\nStart computation of refinement_level " << lvl
                  << " and k = " << k_ << std::flush;

        // start timer
        std::chrono::steady_clock::time_point start_time =
            std::chrono::steady_clock::now();

        // get mesh
        sphere.setRefinementLevel(lvl);
        const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();

        // end timer
        std::chrono::steady_clock::time_point end_time =
            std::chrono::steady_clock::now();
        double elapsed_sec =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                  start_time)
                .count() /
            1000.;

        std::cout << " -> Built Mesh " << elapsed_sec << " [s] " << std::flush;

        solutions_lap.mesh[lvl] = mesh;
        solutions_dirac.mesh[lvl] = mesh;

        // setup the problem
        hodge_builder.SetMesh(mesh);
        hodge_builder.SetK(k_);

        dirac_builder.SetMesh(mesh);
        dirac_builder.SetK(k_);

        // start timer
        start_time = std::chrono::steady_clock::now();

        // compute the system
        hodge_builder.Compute();
        dirac_builder.Compute();

        // end timer
        end_time = std::chrono::steady_clock::now();
        elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_time - start_time)
                          .count() /
                      1000.;

        std::cout << " -> Computed LSE " << elapsed_sec << " [s] "
                  << std::flush;

        // start timer
        start_time = std::chrono::steady_clock::now();

        // solve the system
        hodge_builder.Solve();
        dirac_builder.Solve();

        // end timer
        end_time = std::chrono::steady_clock::now();

        elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_time - start_time)
                          .count() /
                      1000.;

        std::cout << " -> Solved System " << elapsed_sec << " [s] "
                  << std::flush;

        // store solutions
        solutions_lap.solutions[lvl].mu_zero[ik] = hodge_builder.GetMuZero();
        solutions_lap.solutions[lvl].mu_one[ik] =
            std::get<0>(hodge_builder.GetMuOne());
        solutions_lap.solutions[lvl].mu_two[ik] =
            std::get<1>(hodge_builder.GetMuTwo());

        solutions_dirac.solutions[lvl].mu_zero[ik] = dirac_builder.GetMu(0);
        solutions_dirac.solutions[lvl].mu_one[ik] = dirac_builder.GetMu(1);
        solutions_dirac.solutions[lvl].mu_two[ik] = dirac_builder.GetMu(2);

      }  // end loop k
    }    // end loop level

    // end timer total
    std::chrono::steady_clock::time_point end_time_total =
        std::chrono::steady_clock::now();
    double elapsed_sec_total =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_time_total -
                                                              start_time_total)
            .count() /
        1000.;

    std::cout << "\nTotal computation time for all levels " << elapsed_sec_total
              << " [s]\n";

    projects::hldo_sphere::post_processing::compare_results<complex>(
        name_, solutions_lap, solutions_dirac);
  }

 private:
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_zero_;
  std::function<Eigen::Matrix<complex, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f_one_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_two_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_zero_dirac_;
  std::function<Eigen::Matrix<complex, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f_one_dirac_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f_two_dirac_;
  double &k_;
  std::string name_;
};

}  // namespace projects::hldo_sphere::experiments

#endif  // THESIS_EXPERIMENTS_HODGE_AND_DIRAC_H
