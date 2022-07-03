#ifndef HLDO_SPHERE_DEBUGGING_WHITNEY_ONE_BASIS_EXPANSION_COEFFS_H
#define HLDO_SPHERE_DEBUGGING_WHITNEY_ONE_BASIS_EXPANSION_COEFFS_H

/**
 * @file whitney_one_basis_expansion_coeffs.h
 */

#include <lf/assemble/coomatrix.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <norms.h>
#include <results_processing.h>
#include <sphere_triag_mesh_builder.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "debugging.h"

namespace projects::hldo_sphere {

namespace debugging {

using lf::uscalfe::operator-;

/**
 * @brief Class to find the "best" basisexpansion coefficients with the
 * analytical solution
 *
 * The best basis expansion coefficients are defined as the coefficients
 * with the smallest sum of residuals when compared to the analytical solution
 * evaluated at the edge midpoints.
 *
 * The exact  evaluation is not possible because we have three times as many
 * constraints as coefficients when choosing to evaluate at the edge midpoints.
 *
 */
class WhitneyOneBasisExpansionCoeffs {
 public:
  /**
   * @brief Constructor
   *
   * Assigns zero function as analytical solution can be set with set function.
   *
   */
  WhitneyOneBasisExpansionCoeffs() {
    auto u = [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<double, 3, 1> {
      return Eigen::Matrix<double, 3, 1>::Zero();
    };
    u_ = u;
    mu_ = Eigen::VectorXd::Zero(1);
    res_ = Eigen::VectorXd::Zero(1);

    // create mesh factory for basic mesh only because empty values do not
    // complile
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);
    std::shared_ptr<projects::hldo_sphere::mesh::SphereTriagMeshBuilder>
        sphere = std::make_shared<
            projects::hldo_sphere::mesh::SphereTriagMeshBuilder>(
            std::move(factory));
    sphere->setRefinementLevel(0);
    sphere->setRadius(1);
    mesh_ = sphere->Build();
  }

  /**
   * @brief Computes the "best" basisexpansion coefficients.
   *
   * Best in the sense of smallest sum of residuals error when evaluating at the
   * edge midpoints
   *
   * requires the mesh and the analytical functor to be set beforehand.
   */
  void Compute() {
    const lf::assemble::size_type n_edge(mesh_->NumEntities(1));
    lf::assemble::COOMatrix<double> coo_mat(n_edge * 3, n_edge);

    // Compute the Constraint Matrix
    // By iterating over all the triangles and adding the necessary entries
    //
    // That is for each triangle compute the local basis functions and add their
    // contributions at all the edge midpoints to the matrix
    for (const lf::mesh::Entity *cell_p : mesh_->Entities(0)) {
      // Only triangles are supported
      LF_VERIFY_MSG(cell_p->RefEl() == lf::base::RefEl::kTria(),
                    "Unsupported cell type " << cell_p->RefEl());

      // Get the geometry of the entity
      const auto *geom = cell_p->Geometry();

      // Compute the global vertex coordinates
      Eigen::MatrixXd vertices = geom->Global(cell_p->RefEl().NodeCoords());

      // Construct the basis functions from curl of the standard hat functions
      const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
      // The gradients are constant on the triangle
      const Eigen::MatrixXd ref_grads =
          hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
              .transpose();
      // The JacobianInverseGramian is constant on the triangle
      const Eigen::MatrixXd J_inv_trans =
          geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));

      // get the gradients
      const Eigen::MatrixXd grad = J_inv_trans * ref_grads;

      // get edge orientations
      auto edgeOrientations = cell_p->RelativeOrientations();
      Eigen::Vector3d s;
      s << lf::mesh::to_sign(edgeOrientations[0]),
          lf::mesh::to_sign(edgeOrientations[1]),
          lf::mesh::to_sign(edgeOrientations[2]);

      auto edges = cell_p->SubEntities(1);

      Eigen::Matrix<lf::base::size_type, 3, 1> global_idx;
      global_idx << mesh_->Index(*edges[0]), mesh_->Index(*edges[1]),
          mesh_->Index(*edges[2]);

      // loop over edges in the cell
      for (int i = 0; i < 3; i++) {
        int i_p1 = (i + 1) % 3;
        int i_p2 = (i + 2) % 3;

        // compute entries (note we use weight 1/2
        Eigen::Vector3d ent_i_p1 = s(i_p1) / 4. * grad.col(i_p2);
        Eigen::Vector3d ent_i_p2 = -s(i_p2) / 4. * grad.col(i_p2);
        Eigen::Vector3d ent_i = s(i) / 4. * (grad.col(i_p1) - grad.col(i));

        // Add Matrix entries
        coo_mat.AddToEntry(3 * global_idx[i], global_idx[i_p1], ent_i_p1(0));
        coo_mat.AddToEntry(3 * global_idx[i] + 1, global_idx[i_p1],
                           ent_i_p1(1));
        coo_mat.AddToEntry(3 * global_idx[i] + 2, global_idx[i_p1],
                           ent_i_p1(2));

        coo_mat.AddToEntry(3 * global_idx[i], global_idx[i_p2], ent_i_p2(0));
        coo_mat.AddToEntry(3 * global_idx[i] + 1, global_idx[i_p2],
                           ent_i_p2(1));
        coo_mat.AddToEntry(3 * global_idx[i] + 2, global_idx[i_p2],
                           ent_i_p2(2));

        coo_mat.AddToEntry(3 * global_idx[i], global_idx[i], ent_i(0));
        coo_mat.AddToEntry(3 * global_idx[i] + 1, global_idx[i], ent_i(1));
        coo_mat.AddToEntry(3 * global_idx[i] + 2, global_idx[i], ent_i(2));
      }  // loop over local edges
    }    // loop cells

    // compute righthandside vector
    Eigen::VectorXd phi(3 * n_edge);
    for (const lf::mesh::Entity *edge_p : mesh_->Entities(1)) {
      const auto *geom = edge_p->Geometry();
      Eigen::MatrixXd vertices = geom->Global(edge_p->RefEl().NodeCoords());
      Eigen::Vector3d midpoint = 1. / 2. * vertices.col(0) + vertices.col(1);

      // Assign the 3 corresponding entries
      lf::base::size_type glob_idx = mesh_->Index(*edge_p);
      phi.segment(3 * glob_idx, 3) = u_(midpoint);
    }

    // solve the matrix system and store the result
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
        solver;
    Eigen::SparseMatrix<double> sparse_mat = coo_mat.makeSparse();
    sparse_mat.makeCompressed();
    solver.analyzePattern(sparse_mat);
    solver.factorize(sparse_mat);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Could not decompose the matrix");
    }

    mu_ = solver.solve(phi);

    Eigen::VectorXd phi_hat = sparse_mat * mu_;

    res_ = phi_hat - phi;

  }  // Compute

  /**
   * @brief Sets the analytical solution
   * @param u analytical solution vector field
   */
  void SetFunction(std::function<Eigen::Matrix<double, 3, 1>(
                       const Eigen::Matrix<double, 3, 1> &)>
                       u) {
    u_ = u;
  }

  /*
   * @ brief set the mesh for the computations
   *
   * @param mesh pointer to the mesh
   *
   */
  void SetMesh(std::shared_ptr<const lf::mesh::Mesh> mesh) { mesh_ = mesh; }

  /**
   * @brief returns basis expansion coefficients computed in Compute
   */
  Eigen::VectorXd GetMu() { return mu_; }

  /**
   * @brief returns L2 Error between the approximation and the analytical
   * solution
   */
  double GetL2Error() {
    auto square_vector =
        [](Eigen::Matrix<double, Eigen::Dynamic, 1> a) -> double {
      return a.squaredNorm();
    };
    // define quadrule for norms
    lf::quad::QuadRule qr =
        lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 4);
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyOne mf_uh(mu_,
                                                                         mesh_);
    lf::mesh::utils::MeshFunctionGlobal<decltype(u_)> mf_u(u_);
    auto mf_diff = mf_uh - mf_u;
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> l2 =
        projects::hldo_sphere::post_processing::L2norm(mesh_, mf_diff,
                                                       square_vector, qr);
    double error = std::get<0>(l2);
    return error;
  }

  /**
   * @brief returns the squared sum of all residuals at the edge midpoints, this
   * is the quanity minimized in the computation.
   */
  double GetSquaredResidualError() { return res_.squaredNorm(); }

  /*
   * @brief conducts convergence experiment with the given number of
  refinement_level and
   * the function u
   *
   * @param refinement_level the maximal refinement_level considered
   * @param u function to approximate
   *
   */
  static void Experiemnt(
      unsigned refinement_level,
      std::function<Eigen::Vector3d(const Eigen::Vector3d &)> u) {
    std::vector<std::shared_ptr<const lf::mesh::Mesh>> meshs =
        std::vector<std::shared_ptr<const lf::mesh::Mesh>>(refinement_level +
                                                           1);

    // create vecotor with the l2 error in the first component and the
    // squaredresiduals in the second
    std::vector<std::pair<double, double>> solutions(refinement_level + 1);

    // Read the mesh
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

    projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
        projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

    // we always take the same radius
    sphere.setRadius(1.);

    // create a coefficients creator object
    WhitneyOneBasisExpansionCoeffs coeff_creator;
    coeff_creator.SetFunction(u);

    // start timer for total time
    std::chrono::steady_clock::time_point start_time_total =
        std::chrono::steady_clock::now();

    for (unsigned il = 0; il <= refinement_level; ++il) {
      std::cout << "\nStart computation of refinement_level " << il
                << std::flush;

      // start timer
      std::chrono::steady_clock::time_point start_time =
          std::chrono::steady_clock::now();

      // get mesh
      sphere.setRefinementLevel(il);
      const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();
      meshs[il] = mesh;

      coeff_creator.SetMesh(mesh);

      // end timer
      std::chrono::steady_clock::time_point end_time =
          std::chrono::steady_clock::now();
      double elapsed_sec =
          std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                start_time)
              .count() /
          1000.;

      std::cout << " -> Built Mesh " << elapsed_sec << " [s] " << std::flush;

      // start timer
      start_time = std::chrono::steady_clock::now();

      coeff_creator.Compute();

      // end timer
      end_time = std::chrono::steady_clock::now();
      elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time)
                        .count() /
                    1000.;

      std::cout << " -> Computed basis expansion coefficients " << elapsed_sec
                << " [s] " << std::flush;

      // start timer
      start_time = std::chrono::steady_clock::now();

      solutions[il].first = coeff_creator.GetL2Error();
      solutions[il].second = coeff_creator.GetSquaredResidualError();

      // end timer
      end_time = std::chrono::steady_clock::now();

      elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(
                        end_time - start_time)
                        .count() /
                    1000.;

      std::cout << " -> Computed Errors " << elapsed_sec << " [s] "
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

    // create results direcotry
    std::string main_dir = "results";
    std::filesystem::create_directory(main_dir);

    // prepare output
    int table_width = 13;
    int precision = 4;

    // create csv file
    std::ofstream csv_file;
    std::string csv_name = concat("basis_expansion_coeffs_", refinement_level);
    std::replace(csv_name.begin(), csv_name.end(), '.', '_');
    csv_file.open(concat(main_dir, "/", csv_name, ".csv"));
    csv_file << "numCells,"
             << "numEdges,"
             << "numVerts,"
             << "hMax";

    // introduce columns for errors
    csv_file << ",L2Error,SquaredRes";
    csv_file << std::endl;

    Eigen::VectorXd hMax(refinement_level + 1);
    hMax.setZero();

    // loop over all levels contained in the solution
    for (lf::base::size_type lvl = 0; lvl <= refinement_level; ++lvl) {
      // create mesh Functions for solutions the level
      auto &sol_mesh = meshs[lvl];

      // compute meshwidth
      for (const lf::mesh::Entity *e : sol_mesh->Entities(1)) {
        double h = lf::geometry::Volume(*(e->Geometry()));
        if (h > hMax(lvl)) hMax(lvl) = h;
      }

      // print level and mesh informations
      csv_file << sol_mesh->NumEntities(0) << "," << sol_mesh->NumEntities(1)
               << "," << sol_mesh->NumEntities(2) << "," << hMax(lvl);

      csv_file << "," << solutions[lvl].first << "," << solutions[lvl].second
               << "\n";
    }  // end loop over levels

    // close csv file
    csv_file.close();
  }

 private:
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      u_;
  Eigen::VectorXd mu_;
  Eigen::VectorXd res_;
  std::shared_ptr<const lf::mesh::Mesh> mesh_;
};

}  // namespace debugging

}  // namespace projects::hldo_sphere

#endif  // THESIS_DEBUGGING_WHITNEY_ONE_BASIS_EXPANSION_COEFFS_H
