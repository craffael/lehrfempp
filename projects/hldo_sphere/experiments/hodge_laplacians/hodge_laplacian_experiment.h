#ifndef THESIS_EXPERIMENTS_HODGE_LAPLACIAN_H
#define THESIS_EXPERIMENTS_HODGE_LAPLACIAN_H

/**
 * @file hodge_laplacian_experiment.h
 *
 * @brief provides a function for generating solutions, given the validation and
 * load functions
 */
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

/**
 * @brief Creates and solves the discretised Hodge Laplacian source problems for
 * a given list of levels and values of k
 *
 */
class HodgeLaplacianExperiment {
 public:
  /**
   *
   * @brief Constructor setting all the functions and the reference k
   *
   * @param u_zero analytical soltions corresponding to the load function f_zero
   * @param u_one analytical soltions corresponding to the load function f_one
   * @param u_two analytical soltions corresponding to the load function f_two
   * @param f_zero load function corresponding to the analytical solution u_zero
   * @param f_one load function corresponding to the analytical solution u_one
   * @param f_two load function corresponding to the analytical solution u_two
   * @param k reference used in all the functions such that changes of k affect
   * the functions
   * @param name identifier of the example (cretes a folder with this name for
   * the results)
   *
   */
  HodgeLaplacianExperiment(
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> u_zero,
      std::function<
          Eigen::Matrix<double, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          u_one,
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> u_two,
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> f_zero,
      std::function<
          Eigen::Matrix<double, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f_one,
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> f_two,
      double &k, std::string name)
      : u_zero_(u_zero),
        u_one_(u_one),
        u_two_(u_two),
        f_zero_(f_zero),
        f_one_(f_one),
        f_two_(f_two),
        k_(k),
        name_(name) {}

  /**
   *
   * @brief Solves the hodge laplacian source problems for the tensor product of
   * passed refinement levels and ks
   *
   * @param refinement_levels integer list containig all the levels
   * @param ks list of all ks to be used
   *
   */
  void Compute(std::vector<unsigned> refinement_levels,
               std::vector<double> ks) {
    int nl = refinement_levels.size();
    int nk = ks.size();

    // Initialize solution wrapper
    projects::hldo_sphere::post_processing::ProblemSolutionWrapper<double>
        solutions;
    solutions.k = ks;
    solutions.levels = refinement_levels;
    solutions.mesh = std::vector<std::shared_ptr<const lf::mesh::Mesh>>(nl);
    solutions.solutions = std::vector<
        projects::hldo_sphere::post_processing::ProblemSolution<double>>(nl);

    projects::hldo_sphere::operators::HodgeLaplaciansSourceProblems lse_builder;

    // functions are passed by reference hence changing the k still influences
    // the functions
    lse_builder.SetLoadFunctions(f_zero_, f_one_, f_two_);

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
      solutions.solutions[il].mu_zero.resize(nk);
      solutions.solutions[il].mu_one.resize(nk);
      solutions.solutions[il].mu_two.resize(nk);

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

        solutions.mesh[lvl] = mesh;

        // setup the problem
        lse_builder.SetMesh(mesh);
        lse_builder.SetK(k_);

        // start timer
        start_time = std::chrono::steady_clock::now();

        // compute the system
        lse_builder.Compute();

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
        lse_builder.Solve();

        // end timer
        end_time = std::chrono::steady_clock::now();

        elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_time - start_time)
                          .count() /
                      1000.;

        std::cout << " -> Solved System " << elapsed_sec << " [s] "
                  << std::flush;

        // store solutions
        solutions.solutions[lvl].mu_zero[ik] = lse_builder.GetMuZero();
        solutions.solutions[lvl].mu_one[ik] =
            std::get<0>(lse_builder.GetMuOne());
        solutions.solutions[lvl].mu_two[ik] =
            std::get<1>(lse_builder.GetMuTwo());
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

    projects::hldo_sphere::post_processing::process_results<
        decltype(u_zero_), decltype(u_one_), decltype(u_two_), double>(
        name_, solutions, u_zero_, u_one_, u_two_, k_);
  }

 private:
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> u_zero_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      u_one_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> u_two_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f_zero_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f_one_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f_two_;
  double &k_;
  std::string name_;
};

}  // namespace projects::hldo_sphere::experiments

#endif  // THESIS_EXPERIMENTS_HODGE_LAPLACIAN_H
