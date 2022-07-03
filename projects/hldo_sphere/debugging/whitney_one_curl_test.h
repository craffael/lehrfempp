#ifndef HLDO_SPHERE_WHITNEY_ONE_CURL_TEST_H
#define HLDO_SPHERE_WHITNEY_ONE_CURL_TEST_H

/**
 * @file whitney_one_curl_test.h
 */

#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/quad/quad.h>
#include <rot_whitney_one_div_matrix_provider.h>
#include <sphere_triag_mesh_builder.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "whitney_one_basis_expansion_coeffs.h"

namespace projects::hldo_sphere {

namespace debugging {

/**
 * @brief Tests convergence of the below bilinear form on the Mesh to the
 * analytical solution
 *
 * @f[
 *       b(\mathbf{u},v) = \int\limits_{\partial \mathbf{S}} \text{div}\,
 * \mathbf{u} \, v \ dS
 *   @f]
 *
 * For this we compute the Galerkin Matrix B and the values at the basis nodes
 * for the analytical solution and then we expect the following quantity to
 * converge to 0 on mesh refinement
 *
 * @f[
 *  | b(u,v) - u^T \, B \, v | \to 0
 * @f]
 *
 */
class WhitneyOneCurlTest {
 public:
  /**
   * @brief Constructor creates an object with zero testfunctions
   */
  WhitneyOneCurlTest() {
    // create basic function everyting 0 valued by default
    auto v = [](Eigen::Matrix<double, 3, 1> x) -> double { return 0; };
    auto u = [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<double, 3, 1> {
      return Eigen::Matrix<double, 3, 1>::Zero();
    };
    v_ = v;
    u_ = u;
    discrete_sols_ = Eigen::VectorXd::Zero(1);
    meshs_ = std::vector<std::shared_ptr<const lf::mesh::Mesh>>(1);
  }

  /**
   * @brief Computes the error up to the refinemet level 'max_ref'
   *
   * @param max_ref maximal refinement level for which the error is computed
   *
   * The results are then stored in a csv file under
   * resluts/whitney_one_curl_test_`max_ref`.csv
   *
   */
  void Compute(int max_ref) {
    // create matrix provider for the matrix
    projects::hldo_sphere::assemble::RotWhitneyOneDivMatrixProvider
        matrix_div_provider;

    // Read the mesh from the gmsh file
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

    projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
        projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

    // we always take the same radius
    sphere.setRadius(1.);

    discrete_sols_.resize(max_ref + 1);

    meshs_.resize(max_ref + 1);

    // start timer for total time
    std::chrono::steady_clock::time_point start_time_total =
        std::chrono::steady_clock::now();

    // Loop over refinement levels
    for (unsigned lvl = 0; lvl <= max_ref; ++lvl) {
      std::cout << "\nstart computation of max_ref " << lvl << std::flush;

      // start timer
      std::chrono::steady_clock::time_point start_time =
          std::chrono::steady_clock::now();

      // get mesh
      sphere.setrefinementlevel(lvl);
      const std::shared_ptr<lf::mesh::mesh> mesh = sphere.build();
      meshs_[lvl] = mesh;

      // end timer
      std::chrono::steady_clock::time_point end_time =
          std::chrono::steady_clock::now();
      double elapsed_sec =
          std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                start_time)
              .count() /
          1000.;

      std::cout << " -> built mesh " << elapsed_sec << " [s] " << std::flush;

      // start timer
      start_time = std::chrono::steady_clock::now();

      // computes the assemby matrix
      const lf::assemble::dofhandler &dof_handler_edge =
          lf::assemble::uniformfedofhandler(mesh,
                                            {{lf::base::refel::ksegment(), 1}});
      const lf::assemble::dofhandler &dof_handler_cell =
          lf::assemble::uniformfedofhandler(mesh,
                                            {{lf::base::refel::ktria(), 1}});

      // create coo matrix
      const lf::assemble::size_type n_dofs_edge(dof_handler_edge.numdofs());
      const lf::assemble::size_type n_dofs_cell(dof_handler_cell.numdofs());

      // note that the rotwoneformdivelementmatrixprovider returns the
      lf::assemble::coomatrix<double> coo_mat(n_dofs_edge, n_dofs_cell);
      coo_mat.setzero();

      // create matrix with n_dofs_edge_rows and n_dofs_cell columns
      lf::assemble::assemblematrixlocally<lf::assemble::coomatrix<double>>(
          0, dof_handler_cell, dof_handler_edge, matrix_div_provider, coo_mat);

      // compute basis expansion coefficients
      whitneyonebasisexpansioncoeffs one_coeffs;
      one_coeffs.setmesh(mesh);
      one_coeffs.setfunction(u_);
      one_coeffs.compute();
      eigen::matrix<double, eigen::dynamic, 1> u_h(n_dofs_edge);
      u_h = one_coeffs.getmu();

      // get analytical values at the cell midpoints for the expansion
      // coefficients of v
      eigen::matrix<double, eigen::dynamic, 1> v_h(n_dofs_cell);
      for (const lf::mesh::entity *c : mesh->entities(0)) {
        eigen::matrixxd points(3, 3);
        points = lf::geometry::corners(*(c->geometry()));
        eigen::vectorxd midpoint =
            (points.col(0) + points.col(1) + points.col(2)) / 3.;

        // scale up on sphere
        midpoint = midpoint / midpoint.norm();

        unsigned glob_idx = mesh->index(*c);
        v_h(glob_idx) = v_(midpoint);
      }

      discrete_sols_[lvl] = u_h.transpose() * coo_mat.makesparse() * v_h;

      // end timer
      end_time = std::chrono::steady_clock::now();
      elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time)
                        .count() /
                    1000.;

      std::cout << " -> computed discret solution " << elapsed_sec << " [s] "
                << std::flush;

    }  // end loop level

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

    // output the erros in a file

    // create results direcotry
    std::string main_dir = "results";
    std::filesystem::create_directory(main_dir);

    // prepare output
    int table_width = 13;
    int precision = 4;

    // create csv file
    std::ofstream csv_file;
    std::string csv_name = concat("whitney_one_curl_test_", max_ref);
    std::replace(csv_name.begin(), csv_name.end(), '.', '_');
    csv_file.open(concat(main_dir, "/", csv_name, ".csv"));
    csv_file << "numCells,"
             << "numEdges,"
             << "numVerts,"
             << "hMax";

    // introduce columns for errors
    csv_file << ",Errors";
    csv_file << std::endl;

    Eigen::VectorXd hMax(max_ref + 1);
    hMax.setZero();

    // loop over all levels contained in the solution
    for (lf::base::size_type lvl = 0; lvl <= max_ref; ++lvl) {
      // create mesh Functions for solutions the level
      auto &sol_mesh = meshs_[lvl];

      // compute meshwidth
      for (const lf::mesh::Entity *e : sol_mesh->Entities(1)) {
        double h = lf::geometry::Volume(*(e->Geometry()));
        if (h > hMax(lvl)) hMax(lvl) = h;
      }

      // print level and mesh informations
      csv_file << sol_mesh->NumEntities(0) << "," << sol_mesh->NumEntities(1)
               << "," << sol_mesh->NumEntities(2) << "," << hMax(lvl);

      csv_file << "," << std::abs(ana_sol_ - discrete_sols_[lvl]) << "\n";
    }  // end loop over levels

    // close csv file
    csv_file.close();
  }

  /**
   * @brief Sets the value of the analytical solution @f$b(u,v)@f$
   * @param sol value of the analytical soltion
   * @ Note the solution has to comply with the functions u and v, which are
   * passed as functors
   */
  void SetAnaSol(double sol) { ana_sol_ = sol; }

  /**
   * @brief Sets the test functions
   * @param v load function in @f$ L^2 @f$
   * @param u load functions in @f$ H(\text{curl}_{\Gamma}, \partial
   * \mathbf{S})@f$ vector valued
   */
  void SetFunctions(
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> v,
      std::function<
          Eigen::Matrix<double, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          u) {
    v_ = v;
    u_ = u;
  }

 private:
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> v_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      u_;
  double ana_sol_;
  Eigen::VectorXd discrete_sols_;
  std::vector<std::shared_ptr<const lf::mesh::Mesh>> meshs_;
};

}  // namespace debugging

}  // namespace projects::hldo_sphere

#endif  // THESIS_OPERATORS_WHITNEY_ONE_HODGE_LAPLACE_H
