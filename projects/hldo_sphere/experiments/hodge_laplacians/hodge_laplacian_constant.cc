/**
 * @file hodge_laplacian_constant.cc
 * @brief Solve hodge laplacian source problem for fixed function u = 5.2 for
 * the zero and two form  for the zero and two form
 *
 */

#include <hodge_laplacian_experiment.h>

/**
 * @brief Prints the L2 norm errors and the Supremum of the experiment and
 * creates vtk plots for a constant function u = 5.2 for the zero and the two
 * form
 *
 * For the one form we take the tangential vector field @f$ \vec{u}(\vec{x}) =
 * \begin{pmatrix} -x_2 \\ x_1 \\ 0 \end{pmatrix} @f$
 *
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

  // mathematica function output requries the following helpers auto Power =
  // [](double a, double b) -> double { return std::pow(a, b); };
  auto Sin = [](double a) -> double { return std::sin(a); };
  auto Cos = [](double a) -> double { return std::cos(a); };
  auto Sqrt = [](double a) -> double { return std::sqrt(a); };
  auto Power = [](double a, double b) -> double { return std::pow(a, b); };

  double r = 1.;

  // righthandside for the zero and two form
  auto f_zero = [&](const Eigen::Vector3d &x_vec) -> double {
    return k * k * 5.2;
  };

  // righthandside for the one form
  auto f_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::VectorXd {
    // first scale to the sphere
    Eigen::Vector3d x_ = x_vec;
    double x = x_(0);
    double y = x_(1);
    double z = x_(2);

    Eigen::VectorXd ret(3);

    ret << -(Power(k, 2) * y) -
               (2 * y) / (Power(x, 2) + Power(y, 2) + Power(z, 2)),
        x * (Power(k, 2) + 2 / (Power(x, 2) + Power(y, 2) + Power(z, 2))), 0;
    return ret;
  };

  // Compute the analytic solution of the problem
  auto u_zero = [&](const Eigen::Vector3d x_vec) -> double { return 5.2; };

  // Compute the analytic solution of the problem
  auto u_one = [&](const Eigen::Vector3d x_vec) -> Eigen::Vector3d {
    // first scale to the sphere
    Eigen::Vector3d x_ = x_vec;
    double x = x_(0);
    double y = x_(1);
    double z = x_(2);

    Eigen::VectorXd ret(3);

    ret << -y, x, 0;

    return ret;
  };

  projects::hldo_sphere::experiments::HodgeLaplacianExperiment experiment(
      u_zero, u_one, u_zero, f_zero, f_one, f_zero, k, "constant");

  experiment.Compute(refinement_levels, ks);

  return 0;
}
