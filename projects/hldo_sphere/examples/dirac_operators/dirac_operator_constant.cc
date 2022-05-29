/**
 * @file hodge_laplacian_constant.cc
 * @brief Solve hodge laplacian source problem for fixed function u = 5.2
 *
 */

#include <dirac_operator_example.h>

using lf::uscalfe::operator-;
using complex = std::complex<double>;

/**
 * @brief Prints the L2 norm errors and the Supremum of the experiment and
 * creates vtk plots for a constant function u = 5.2
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

  projects::hldo_sphere::examples::DiracOperatorExample example(
      u_zero, u_one, u_two, f_zero, f_one, f_two, k, "constant");

  example.Compute(refinement_levels, ks);

  return 0;
}
