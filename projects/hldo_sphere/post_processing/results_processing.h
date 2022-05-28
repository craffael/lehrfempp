#ifndef THESIS_POST_PROCESSING_RESULTS_PROCESSING_H
#define THESIS_POST_PROCESSING_RESULTS_PROCESSING_H

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
#include <two_form_to_mesh_data_set.h>

#include <Eigen/Dense>
#include <cstdlib>
#include <cstring>
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
 * @brief stores information to recover convergence properties
 *
 * @tparam U_ZERO functor for the analytical solution of the zero form
 * @tparam U_ONE functor for the analytical solution of the one form
 * @tparam U_TWO functor for the analytical solution of the two form
 */
template <typename SCALAR>
struct ProblemSolution {
  std::shared_ptr<const lf::mesh::Mesh> mesh;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> mu_zero;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> mu_one;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> mu_two;
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
 * @brief Generate vtk files containing meshdata, plots and tables with results
 *
 * @tparam U_ZERO functor for the analytical solution of the zero form
 * @tparam U_ONE functor for the analytical solution of the one form
 * @tparam U_TWO functor for the analytical solution of the two form
 * @tparam SCALAR type of the resultvectors
 *
 * @param name identifying the experiment (used to name the outputs)
 * @param results contains a list of ProblemSolution structs, each containing a
 * mesh with analyitical solutions and the coefficient vectors for the
 * approximate solution
 *
 */
template <typename U_ZERO, typename U_ONE, typename U_TWO, typename SCALAR>
void process_results(std::string name,
                     std::vector<ProblemSolution<SCALAR>> &results,
                     U_ZERO &u_zero, U_ONE u_one, U_TWO u_two) {
  int n = results.size();

  // prepare output
  int table_width = 15;
  int precision = 4;

  // create csv file
  std::ofstream csv_file;
  csv_file.open(concat("result_", name, ".csv"));
  csv_file << "numCells,"
           << "numEdges,"
           << "numVerts,"
           << "hMax,"
           << "SupErrorZero,"
           << "SupErrorOne,"
           << "SupErrorTwo,"
           << "L2ErrorZero,"
           << "L2ErrorOne,"
           << "L2ErrorTwo\n";

  // print all the errors
  std::cout << std::endl
            << std::endl
            << std::endl
            << std::setw(table_width) << "Norm Type" << std::setw(table_width)
            << "Ref. level" << std::setw(table_width) << "Num Cells"
            << std::setw(table_width) << "Error Zero" << std::setw(table_width)
            << "Error One" << std::setw(table_width) << "Error Two" << std::endl
            << std::endl;

  // loop over all levels contained in the solution
  for (lf::base::size_type lvl = 0; lvl < n; ++lvl) {
    // create mesh Functions for solutions
    auto &sol = results[lvl];
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyZero mf_zero(
        sol.mu_zero, sol.mesh);
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyOne mf_one(
        sol.mu_one, sol.mesh);
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyTwo mf_two(
        sol.mu_two, sol.mesh);

    lf::mesh::utils::MeshFunctionGlobal<U_ZERO> mf_zero_ana(u_zero);
    lf::mesh::utils::MeshFunctionGlobal<U_ONE> mf_one_ana(u_one);
    lf::mesh::utils::MeshFunctionGlobal<U_TWO> mf_two_ana(u_two);

    // Compute the error of the solutions
    auto mf_diff_zero = mf_zero - mf_zero_ana;
    auto mf_diff_one = mf_one - mf_one_ana;
    auto mf_diff_two = mf_two - mf_zero_ana;

    // Perform post processing on the data
    lf::io::VtkWriter writer(sol.mesh,
                             concat("result_", name, "_", lvl, ".vtk"));

    // define square functions
    auto square_scalar = [](SCALAR a) -> double {
      return std::abs(a) * std::abs(a);
    };
    auto square_vector =
        [](Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> a) -> double {
      return a.squaredNorm();
    };

    // define quadrule for norms
    lf::quad::QuadRule qr =
        lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2);

    // get error for zero form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_zero =
        projects::hldo_sphere::post_processing::L2norm(sol.mesh, mf_diff_zero,
                                                       square_scalar, qr);
    const double error_zero = std::get<0>(L2_zero);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_zero =
        std::get<1>(L2_zero);
    writer.WriteCellData(concat("mf_zero_diff_", lvl), data_set_error_zero);

    // get error for one form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_one =
        projects::hldo_sphere::post_processing::L2norm(sol.mesh, mf_diff_one,
                                                       square_vector, qr);
    const double error_one = std::get<0>(L2_one);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_one =
        std::get<1>(L2_one);
    writer.WriteCellData(concat("mf_one_diff_", lvl), data_set_error_one);

    // get error for two form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_two =
        projects::hldo_sphere::post_processing::L2norm(sol.mesh, mf_diff_two,
                                                       square_scalar, qr);
    const double error_two = std::get<0>(L2_two);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_two =
        std::get<1>(L2_two);
    writer.WriteCellData(concat("mf_two_diff_", lvl), data_set_error_two);

    // get sup error for zero form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>>
        sup_zero = projects::hldo_sphere::post_processing::SupNorm(
            sol.mesh, mf_diff_zero, square_scalar, qr);
    const double error_sup_zero = std::get<0>(sup_zero);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_sup_zero =
        std::get<1>(sup_zero);
    writer.WriteCellData(concat("mf_zero_diff_sup_", lvl),
                         data_set_error_sup_zero);

    // get error for one form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> sup_one =
        projects::hldo_sphere::post_processing::SupNorm(sol.mesh, mf_diff_one,
                                                        square_vector, qr);
    const double error_sup_one = std::get<0>(sup_one);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_sup_one =
        std::get<1>(sup_one);
    writer.WriteCellData(concat("mf_one_diff_sup_", lvl),
                         data_set_error_sup_one);

    // get error for two form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> sup_two =
        projects::hldo_sphere::post_processing::SupNorm(sol.mesh, mf_diff_two,
                                                        square_scalar, qr);
    const double error_sup_two = std::get<0>(sup_two);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_sup_two =
        std::get<1>(sup_two);
    writer.WriteCellData(concat("mf_two_diff_sup_", lvl),
                         data_set_error_sup_two);

    // print all the L2 errors
    std::cout << std::setprecision(precision) << std::setw(table_width)
              << "L2_norm^2" << std::setw(table_width) << lvl
              << std::setw(table_width) << sol.mesh->NumEntities(0)
              << std::setw(table_width) << error_zero << std::setw(table_width)
              << error_one << std::setw(table_width) << error_two << std::endl;

    // print all the sup errors
    std::cout << std::setw(table_width) << "sup_norm^2"
              << std::setw(table_width) << lvl << std::setw(table_width)
              << sol.mesh->NumEntities(0) << std::setw(table_width)
              << error_sup_zero << std::setw(table_width) << error_sup_one
              << std::setw(table_width) << error_sup_two << std::endl
              << std::endl;

    /******************************
     * write new line in csv file
     ******************************/

    // compute meshwidth
    double h_max = 0;
    for (const lf::mesh::Entity *e : sol.mesh->Entities(1)) {
      double h = lf::geometry::Volume(*(e->Geometry()));
      if (h > h_max) h_max = h;
    }
    csv_file << sol.mesh->NumEntities(0) << "," << sol.mesh->NumEntities(1)
             << "," << sol.mesh->NumEntities(2) << "," << h_max << ","
             << error_sup_zero << "," << error_sup_one << "," << error_sup_two
             << "," << error_zero << "," << error_one << "," << error_two
             << "\n";
  }

  // close csv file
  csv_file.close();
}

}  // end namespace post_processing

}  // namespace projects::hldo_sphere

#endif  // THESIS_POST_PROCESSING_RESULTS_PROCESSING_H
