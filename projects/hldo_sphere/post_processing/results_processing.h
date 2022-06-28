#ifndef HLDO_SPHERE_POST_PROCESSING_RESULTS_PROCESSING_H
#define HLDO_SPHERE_POST_PROCESSING_RESULTS_PROCESSING_H

/**
 * @file result_processing
 *
 * @brief takes a list of results and creates vtk files and tables
 */

#include <lf/io/vtk_writer.h>
#include <mesh_function_whitney_one.h>
#include <mesh_function_whitney_two.h>
#include <mesh_function_whitney_zero.h>
#include <norms.h>
#include <sphere_triag_mesh_builder.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace projects::hldo_sphere {

namespace post_processing {

// operator used to subract meshfunctions
using lf::uscalfe::operator-;

/**
 * @brief stores solutions for several k and one refinement level
 *
 * @tparam SCALAR base type of the solutions vectors
 */
template <typename SCALAR>
struct ProblemSolution {
  std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>> mu_zero;
  std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>> mu_one;
  std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>> mu_two;
};

/**
 * @brief stores solutions and setupinformation for a list of refinement levels
 * and k values
 *
 * @tparam SCALAR base type of the solutions vectors
 */
template <typename SCALAR>
struct ProblemSolutionWrapper {
  std::vector<unsigned> levels;
  std::vector<double> k;
  // each element in the vector corresponds to one refinement level
  std::vector<std::shared_ptr<const lf::mesh::Mesh>> mesh;
  std::vector<ProblemSolution<SCALAR>> solutions;
};

/**
 * @brief Concatenate objects defining an operator<<(std::ostream&)
 * @param args A variadic pack of objects implementing
 * `operator<<(std::ostream&)`
 * @returns A string with the objects concatenated
 */
template <typename... Args>
static std::string concat(Args &&...args) {
  std::ostringstream ss;
  (ss << ... << args);
  return ss.str();
}

/**
 * @brief Generate vtk files containing meshdata and a csv table with results
 *
 * The Table contains for each refinement level, the number of cells, edges and
 * vetices, as well as the global mesh width (the maximum edgelenght of a
 * triangle). Moreover it contains the l2 error and the Supremum error as well
 * as the approximate orders for every value of k and all three forms.
 *
 *
 * @tparam U_ZERO functor for the analytical solution of the zero form
 * @tparam U_ONE functor for the analytical solution of the one form
 * @tparam U_TWO functor for the analytical solution of the two form
 * @tparam SCALAR base type of the resultvectors
 *
 * @param name identifying the experiment (used to name the outputs)
 * @param results ProblemSolutionWrapper with the results for all levels and ks
 * @param u_zero vector with the analytical solutions for u_zero of every k
 * @param u_one vector with the analytical solutions for u_one of every k
 * @param u_two vector with the analytical solutions for u_two of every k
 * @param k reference to the variable used in the functions u_one, u_two, u_zero
 *
 */
template <typename U_ZERO, typename U_ONE, typename U_TWO, typename SCALAR>
void process_results(std::string name, ProblemSolutionWrapper<SCALAR> &results,
                     U_ZERO &u_zero, U_ONE &u_one, U_TWO &u_two, double &k) {
  // create results direcotry
  std::string main_dir = "results";
  std::filesystem::create_directory(main_dir);

  // each containing number of k values
  int nk = results.k.size();
  int nl = results.levels.size();

  // prepare output
  int table_width = 13;
  int precision = 4;

  // name the csv file accoring to min and max k
  double min_k = results.k[0];
  double max_k = results.k[nk - 1];

  // create csv file
  std::ofstream csv_file;
  std::string csv_name = concat("result_", name, "_", min_k, "-", max_k);
  std::replace(csv_name.begin(), csv_name.end(), '.', '_');
  csv_file.open(concat(main_dir, "/", csv_name, ".csv"));
  csv_file << "numCells,"
           << "numEdges,"
           << "numVerts,"
           << "hMax";

  // define square functions for norms
  auto square_scalar = [](SCALAR a) -> double {
    return std::abs(a) * std::abs(a);
  };
  auto square_vector =
      [](Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> a) -> double {
    return a.squaredNorm();
  };

  // define quadrule for norms
  lf::quad::QuadRule qr = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 4);

  // create solution matrix with levels on rows and ks in colums
  Eigen::MatrixXd SupErrorZero = Eigen::MatrixXd::Zero(nl, nk);
  Eigen::MatrixXd SupErrorOne = Eigen::MatrixXd::Zero(nl, nk);
  Eigen::MatrixXd SupErrorTwo = Eigen::MatrixXd::Zero(nl, nk);
  Eigen::MatrixXd L2ErrorZero = Eigen::MatrixXd::Zero(nl, nk);
  Eigen::MatrixXd L2ErrorOne = Eigen::MatrixXd::Zero(nl, nk);
  Eigen::MatrixXd L2ErrorTwo = Eigen::MatrixXd::Zero(nl, nk);
  Eigen::VectorXd hMax = Eigen::VectorXd::Zero(nl);

  // tabulate every k value
  for (int i = 0; i < nk; i++) {
    csv_file << ",SupErrorZero_" << results.k[i] << ",RateSupErrorZero_"
             << results.k[i] << ",SupErrorOne_" << results.k[i]
             << ",RateSupErrorOne_" << results.k[i] << ",SupErrorTwo_"
             << results.k[i] << ",RateSupErrorTwo_" << results.k[i]
             << ",L2ErrorZero_" << results.k[i] << ",RateL2ErrorZero_"
             << results.k[i] << ",L2ErrorOne_" << results.k[i]
             << ",RateL2ErrorOne_" << results.k[i] << ",L2ErrorTwo_"
             << results.k[i] << ",RateL2ErrorTwo_" << results.k[i];

    // create direcotries for each k
    std::string kstr = concat(results.k[i]);
    std::replace(kstr.begin(), kstr.end(), '.', '_');
    std::filesystem::create_directory(concat(main_dir, "/k_", kstr));
  }
  csv_file << std::endl;

  // write the k values to the csv file
  csv_file << ","
           << ","
           << ",";
  for (int i = 0; i < nk; i++) {
    csv_file << "," << results.k[i] << "," << results.k[i] << ","
             << results.k[i] << "," << results.k[i] << "," << results.k[i]
             << "," << results.k[i] << "," << results.k[i] << ","
             << results.k[i] << "," << results.k[i] << "," << results.k[i]
             << "," << results.k[i] << "," << results.k[i];
  }
  csv_file << std::endl;

  // loop over all levels contained in the solution
  for (lf::base::size_type lvl = 0; lvl < nl; ++lvl) {
    // create mesh Functions for solutions the level
    auto &sol_mesh = results.mesh[lvl];
    auto &sol_mus = results.solutions[lvl];

    // compute meshwidth
    for (const lf::mesh::Entity *e : sol_mesh->Entities(1)) {
      double h = lf::geometry::Volume(*(e->Geometry()));
      if (h > hMax(lvl)) hMax(lvl) = h;
    }

    // print level and mesh informations
    csv_file << sol_mesh->NumEntities(0) << "," << sol_mesh->NumEntities(1)
             << "," << sol_mesh->NumEntities(2) << "," << hMax(lvl);

    // compute the errors for each k
    for (int ik = 0; ik < nk; ik++) {
      std::string kstr = concat(results.k[ik]);
      std::replace(kstr.begin(), kstr.end(), '.', '_');
      std::string folder = concat(main_dir, "/k_", kstr);

      projects::hldo_sphere::post_processing::MeshFunctionWhitneyZero mf_zero(
          sol_mus.mu_zero[ik], sol_mesh);
      projects::hldo_sphere::post_processing::MeshFunctionWhitneyOne mf_one(
          sol_mus.mu_one[ik], sol_mesh);
      projects::hldo_sphere::post_processing::MeshFunctionWhitneyTwo mf_two(
          sol_mus.mu_two[ik], sol_mesh);

      // set k and therefore change function since the k is used in the
      // funcitons
      k = results.k[ik];
      lf::mesh::utils::MeshFunctionGlobal<U_ZERO> mf_zero_ana(u_zero);
      lf::mesh::utils::MeshFunctionGlobal<U_ONE> mf_one_ana(u_one);
      lf::mesh::utils::MeshFunctionGlobal<U_TWO> mf_two_ana(u_two);

      // Compute the error of the solutions
      auto mf_diff_zero = mf_zero - mf_zero_ana;
      auto mf_diff_one = mf_one - mf_one_ana;
      auto mf_diff_two = mf_two - mf_two_ana;

      // Perform post processing on the data
      lf::io::VtkWriter writer(sol_mesh, concat(folder, "/result_", name, "_",
                                                lvl, "_k_", kstr, ".vtk"));

      // get error for zero form
      const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
          L2_zero = projects::hldo_sphere::post_processing::L2norm(
              sol_mesh, mf_diff_zero, square_scalar, qr);
      L2ErrorZero(lvl, ik) = std::get<0>(L2_zero);
      const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_zero =
          std::get<1>(L2_zero);
      writer.WriteCellData(concat("mf_zero_diff_", lvl), data_set_error_zero);

      // get error for one form
      const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
          L2_one = projects::hldo_sphere::post_processing::L2norm(
              sol_mesh, mf_diff_one, square_vector, qr);
      L2ErrorOne(lvl, ik) = std::get<0>(L2_one);
      const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_one =
          std::get<1>(L2_one);
      writer.WriteCellData(concat("mf_one_diff_", lvl), data_set_error_one);

      // get error for two form
      const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
          L2_two = projects::hldo_sphere::post_processing::L2norm(
              sol_mesh, mf_diff_two, square_scalar, qr);
      L2ErrorTwo(lvl, ik) = std::get<0>(L2_two);
      const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_two =
          std::get<1>(L2_two);
      writer.WriteCellData(concat("mf_two_diff_", lvl), data_set_error_two);

      // get sup error for zero form
      const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
          sup_zero = projects::hldo_sphere::post_processing::SupNorm(
              sol_mesh, mf_diff_zero, square_scalar, qr);
      SupErrorZero(lvl, ik) = std::get<0>(sup_zero);
      const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_sup_zero =
          std::get<1>(sup_zero);
      writer.WriteCellData(concat("mf_zero_diff_sup_", lvl),
                           data_set_error_sup_zero);

      // get error for one form
      const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
          sup_one = projects::hldo_sphere::post_processing::SupNorm(
              sol_mesh, mf_diff_one, square_vector, qr);
      SupErrorOne(lvl, ik) = std::get<0>(sup_one);
      const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_sup_one =
          std::get<1>(sup_one);
      writer.WriteCellData(concat("mf_one_diff_sup_", lvl),
                           data_set_error_sup_one);

      // get error for two form
      const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
          sup_two = projects::hldo_sphere::post_processing::SupNorm(
              sol_mesh, mf_diff_two, square_scalar, qr);
      SupErrorTwo(lvl, ik) = std::get<0>(sup_two);
      const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_sup_two =
          std::get<1>(sup_two);
      writer.WriteCellData(concat("mf_two_diff_sup_", lvl),
                           data_set_error_sup_two);

      double l2RateZero = 0;
      double l2RateOne = 0;
      double l2RateTwo = 0;
      double supRateZero = 0;
      double supRateOne = 0;
      double supRateTwo = 0;

      // we can only compute order form the second level
      if (lvl > 0) {
        l2RateZero =
            (log(L2ErrorZero(lvl, ik)) - log(L2ErrorZero(lvl - 1, ik))) /
            (log(hMax(lvl) / hMax(lvl - 1)));
        l2RateOne = (log(L2ErrorOne(lvl, ik)) - log(L2ErrorOne(lvl - 1, ik))) /
                    (log(hMax(lvl) / hMax(lvl - 1)));
        l2RateTwo = (log(L2ErrorTwo(lvl, ik)) - log(L2ErrorTwo(lvl - 1, ik))) /
                    (log(hMax(lvl) / hMax(lvl - 1)));
        supRateZero =
            (log(SupErrorZero(lvl, ik)) - log(SupErrorZero(lvl - 1, ik))) /
            (log(hMax(lvl) / hMax(lvl - 1)));
        supRateOne =
            (log(SupErrorOne(lvl, ik)) - log(SupErrorOne(lvl - 1, ik))) /
            (log(hMax(lvl) / hMax(lvl - 1)));
        supRateTwo =
            (log(SupErrorTwo(lvl, ik)) - log(SupErrorTwo(lvl - 1, ik))) /
            (log(hMax(lvl) / hMax(lvl - 1)));
      }

      /******************************
       * append solution of the current k in the outputs
       ******************************/

      csv_file << "," << SupErrorZero(lvl, ik) << "," << supRateZero << ","
               << SupErrorOne(lvl, ik) << "," << supRateOne << ","
               << SupErrorTwo(lvl, ik) << "," << supRateTwo << ","
               << L2ErrorZero(lvl, ik) << "," << l2RateZero << ","
               << L2ErrorOne(lvl, ik) << "," << l2RateOne << ","
               << L2ErrorTwo(lvl, ik) << "," << l2RateTwo;

    }  // end loop over k

    csv_file << "\n";
  }  // end loop over levels

  // close csv file
  csv_file.close();
}

}  // end namespace post_processing

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_POST_PROCESSING_RESULTS_PROCESSING_H
