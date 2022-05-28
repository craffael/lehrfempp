/**
 * @file hodge_laplacian_constant.cc
 * @brief Solve hodge laplacian source problem for fixed function u = 5.2
 *
 */

#include <dirac_operator_source_problem.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/utils.h>
#include <norms.h>
#include <results_processing.h>
#include <sphere_triag_mesh_builder.h>

#include <chrono>

using lf::uscalfe::operator-;
using complex = std::complex<double>;

/**
 * @brief Prints the L2 norm errors and the Supremum of the experiment and
 * creates vtk plots for a constant function u = 5.2
 */
int main(int argc, char *argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] << " max_refinement_level k radius"
              << std::endl;
    exit(1);
  }

  double k;
  const unsigned refinement_level = atoi(argv[1]);
  double r;
  sscanf(argv[2], "%lf", &k);
  sscanf(argv[3], "%lf", &r);
  std::cout << "max_refinement_level : " << refinement_level << std::endl;

  // Read the mesh from the gmsh file
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  // we always take the same radius
  sphere.setRadius(r);

  // righthandside for the zero
  auto f_zero = [&](const Eigen::Vector3d &x_vec) -> complex {
    return complex(0, k * 5.2);
  };

  // righthandside for the one form
  auto f_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::VectorXcd {
    Eigen::Vector3cd ret;
    ret << complex(0, -k * x_vec(1)), complex(0, k * x_vec(0)), 0;
    return ret;
  };

  auto Power = [](double a, double b) -> double { return pow(a, b); };
  auto Sqrt = [](double a) -> double { return sqrt(a); };
  auto Complex = [](double a, double b) -> complex { return complex(a, b); };

  // righthandside for the two form
  auto f_two = [&](const Eigen::Vector3d &x_vec) -> complex {
    double x = x_vec(0);
    double y = x_vec(1);
    double z = x_vec(2);
    return (Power(x, 8) *
                (-2. * z + Complex(0., 5.2) * k *
                               Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(y, 8) *
                (-2. * z + Complex(0., 5.2) * k *
                               Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(z, 8) *
                (-2. * z + Complex(0., 5.2) * k *
                               Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(y, 6) * Power(z, 2) *
                (-8. * z + Complex(0., 20.8) * k *
                               Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(y, 2) * Power(z, 6) *
                (-8. * z + Complex(0., 20.8) * k *
                               Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(x, 6) * (Power(y, 2) + Power(z, 2)) *
                (-8. * z + Complex(0., 20.8) * k *
                               Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(y, 4) * Power(z, 4) *
                (-12. * z + Complex(0., 31.200000000000003) * k *
                                Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
            Power(x, 4) *
                (Power(y, 4) * (-12. * z + Complex(0., 31.200000000000003) * k *
                                               Sqrt(Power(x, 2) + Power(y, 2) +
                                                    Power(z, 2))) +
                 Power(z, 4) * (-12. * z + Complex(0., 31.200000000000003) * k *
                                               Sqrt(Power(x, 2) + Power(y, 2) +
                                                    Power(z, 2))) +
                 Power(y, 2) * Power(z, 2) *
                     (-24. * z +
                      Complex(0., 62.400000000000006) * k *
                          Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2)))) +
            Power(x, 2) *
                (Power(y, 6) * (-8. * z + Complex(0., 20.8) * k *
                                              Sqrt(Power(x, 2) + Power(y, 2) +
                                                   Power(z, 2))) +
                 Power(z, 6) * (-8. * z + Complex(0., 20.8) * k *
                                              Sqrt(Power(x, 2) + Power(y, 2) +
                                                   Power(z, 2))) +
                 Power(y, 4) * Power(z, 2) *
                     (-24. * z +
                      Complex(0., 62.400000000000006) * k *
                          Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))) +
                 Power(y, 2) * Power(z, 4) *
                     (-24. * z +
                      Complex(0., 62.400000000000006) * k *
                          Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2))))) /
           Power(Power(x, 2) + Power(y, 2) + Power(z, 2), 4.5);
  };

  // Compute the analytic solution of the problem
  auto u_zero = [&](const Eigen::Vector3d x_vec) -> complex { return 5.2; };
  auto u_two = [&](const Eigen::Vector3d x_vec) -> complex { return 5.2; };

  // Compute the analytic solution of the problem
  auto u_one = [&](const Eigen::Vector3d x_vec) -> Eigen::Vector3cd {
    double x = x_vec(0);
    double y = x_vec(1);
    double z = x_vec(2);
    Eigen::Vector3cd ret;
    ret << -y, x, 0;
    return ret;
  };

  // Solve the problem for each mesh in the hierarchy
  std::vector<projects::hldo_sphere::post_processing::ProblemSolution<complex>>
      solutions(refinement_level + 1);

  projects::hldo_sphere::discretization::DiracOperatorSourceProblem lse_builder;
  lse_builder.SetLoadFunctions(f_zero, f_one, f_two);

  // start timer for total time
  std::chrono::steady_clock::time_point start_time_total =
      std::chrono::steady_clock::now();

  for (lf::base::size_type lvl = 0; lvl < refinement_level + 1; ++lvl) {
    std::cout << "\nStart computation of refinement_level " << lvl << " "
              << std::flush;

    // start timer
    std::chrono::steady_clock::time_point start_time =
        std::chrono::steady_clock::now();

    // get mesh
    sphere.setRefinementLevel(lvl);
    const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();

    // end timer
    std::chrono::steady_clock::time_point end_time =
        std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_time - start_time)
                             .count() /
                         1000.;

    std::cout << " -> Built Mesh " << elapsed_sec << " [s] " << std::flush;

    auto &sol = solutions[lvl];
    sol.mesh = mesh;

    // setup the problem
    lse_builder.SetMesh(mesh);
    lse_builder.SetK(k);

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

    std::cout << " -> Computed LSE " << elapsed_sec << " [s] " << std::flush;

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

    std::cout << " -> Solved System " << elapsed_sec << " [s] " << std::flush;

    // store solutions
    sol.mu_zero = lse_builder.GetMu(0);
    sol.mu_one = lse_builder.GetMu(1);
    sol.mu_two = lse_builder.GetMu(2);
  }

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
      decltype(u_zero), decltype(u_one), decltype(u_two), complex>(
      "constant", solutions, u_zero, u_one, u_two);

  return 0;
}
