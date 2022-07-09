/**
 * @file hodge_laplacians.cc
 * @brief Solve hodge laplacian source problem for fixed function u = 0 and
 * variable k
 */

#include <hodge_laplacian_experiment.h>

/**
 * @brief Prints the L2 norm errors and the Supremum of the experiment and
 * creates vtk plots for a constant function u = 0
 */
int main(int argc, char *argv[]) {
  if (argc != 4 && argc != 5) {
    std::cerr << "Usage: " << argv[0]
              << " max_refinement_level min_k max_k step=0.1 " << std::endl;
    exit(1);
  }

  double min_k;
  double max_k;
  double step = 0.1;
  const unsigned refinement_level = atoi(argv[1]);
  sscanf(argv[2], "%lf", &min_k);
  sscanf(argv[3], "%lf", &max_k);
  if (argc == 5) sscanf(argv[4], "%lf", &step);
  std::cout << "max_refinement_level : " << refinement_level << std::endl;
  std::cout << "min_k : " << min_k << ", max_k = " << max_k
            << " step = " << step << std::endl;

  // needs to be a reference such that it can be reassigned
  double &k = min_k;

  std::vector<unsigned> refinement_levels(refinement_level + 1);
  for (int i = 0; i < refinement_level + 1; i++) {
    refinement_levels[i] = i;
  }

  std::vector<double> ks;
  for (double i = min_k; i <= max_k; i += step) {
    ks.push_back(i);
  }

  // righthandside for the zero and two form
  auto f_zero = [&](const Eigen::Vector3d &x_vec) -> double { return 0; };

  // righthandside for the one form
  auto f_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::VectorXd {
    Eigen::VectorXd ret(3);
    ret << 0, 0, 0;
    return ret;
  };

  // Compute the analytic solution of the problem
  auto u_zero = [&](const Eigen::Vector3d &x_vec) -> double { return 0; };

  // Compute the analytic solution of the problem
  auto u_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::Vector3d {
    Eigen::VectorXd ret(3);

    ret << 0, 0, 0;
    return ret;
  };

  projects::hldo_sphere::experiments::HodgeLaplacianExperiment<double>
      experiment(u_zero, u_one, u_zero, f_zero, f_one, f_zero, k, "zero");

  experiment.Compute(refinement_levels, ks);

  return 0;
}
