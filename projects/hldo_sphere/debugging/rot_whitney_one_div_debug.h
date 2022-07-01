#ifndef HLDO_SPHERE_ROT_WHITNEY_ONE_DIV_DEBUG_H
#define HLDO_SPHERE_ROT_WHITNEY_ONE_DIV_DEBUG_H

/**
 * @file rot_whitney_one_div_debug.h
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

namespace projects::hldo_sphere {

namespace debugging {

using SCALAR = std::complex<double>;

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
class RotWhitneyOneDivDebug {
 public:
  /**
   * @brief Constructor
   * creates basic mesh (Octaeder with radius 1.0)
   * creates zerovalued function f
   *
   */
  RotWhitneyOneDivDebug() {
    // create basic function everyting 0 valued by default
    auto v = [](Eigen::Matrix<double, 3, 1> x) -> SCALAR { return 0; };
    auto u = [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<SCALAR, 3, 1> {
      return Eigen::Matrix<SCALAR, 3, 1>::Zero();
    };
    v_ = v;
    u_ = u;
  }

  /**
   * @brief Computes the error up to the refinemet level 'max_ref'
   *
   * @param max_ref maximal refinement level for which the error is computed
   *
   */
  void Compute() {
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

    // Loop over refinement levels
    for (unsigned lvl = 0; lvl < max_ref; ++lvl) {
      std::cout << "\nStart computation of refinement_level " << lvl
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
      double elapsed_sec =
          std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                start_time)
              .count() /
          1000.;

      std::cout << " -> Built Mesh " << elapsed_sec << " [s] " << std::flush;

      // start timer
      start_time = std::chrono::steady_clock::now();

      // Computes the assemby matrix
      const lf::assemble::DofHandler &dof_handler_edge =
          lf::assemble::UniformFEDofHandler(mesh,
                                            {{lf::base::RefEl::kSegment(), 1}});
      const lf::assemble::DofHandler &dof_handler_cell =
          lf::assemble::UniformFEDofHandler(mesh,
                                            {{lf::base::RefEl::kTria(), 1}});

      // create COO Matrix
      const lf::assemble::size_type n_dofs_edge(dof_handler_edge.NumDofs());
      const lf::assemble::size_type n_dofs_cell(dof_handler_cell.NumDofs());

      // note that the RotWOneFormDivElementMatrixProvider returns the
      lf::assemble::COOMatrix<double> coo_mat(n_dofs_edge, n_dofs_cell);
      coo_mat.setZero();

      // create matrix with n_dof_edge_rows and n_dof_cell columns
      lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
          0, dof_handler_cell, dof_handler_edge, matrix_div_provider, coo_mat);

      // compute basis expansion coefficients
      Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> u_h(n_dof_edge);
      u_h.setZero();

      // get analytical values at the edge midpoints
      for (const lf::mesh::Entity e : mesh->Entities(1)) {
        Eigen::MatrixXd points(3, 2);
        points = lf::geometr::Corners(e.Geometry());
        Eigen::VectorXd midpoint = (points(0) + points(1)) / 2.;

        // scale up on sphere
        midpoint = midpoint / midpoint.norm();

        // assign theoretical basis expansion coefficient
        // This would require to solve a linear system of equations with more
        // equations than unknowns which does not necessarily have a solution
        // and hence is not so easily possible
      }

      // get the cell midpoints of the mesh

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
      solutions.solutions[lvl].mu_zero[ik] = lse_builder.GetMu(0);
      solutions.solutions[lvl].mu_one[ik] = lse_builder.GetMu(1);
      solutions.solutions[lvl].mu_two[ik] = lse_builder.GetMu(2);
    }  // end loop k
  }    // end loop level
}

  /**
   * @brief Sets the value of the analytical solution @f$b(u,v)@f$
   * @param sol value of the analytical soltion
   * @ Note the solution has to comply with the functions u and v, which are
   * passed as functors
   */
  void SetAnaSol(SCALAR sol) {
  ana_sol_ = sol;
}

/**
 * @brief Sets the test functions
 * @param v load function in @f$ L^2 @f$
 * @param u load functions in @f$ H(\text{curl}_{\Gamma}, \partial
 * \mathbf{S})@f$ vector valued
 */
void SetFunctions(std::function<SCALAR(const Eigen::Matrix<double, 3, 1> &)> v,
                  std::function<Eigen::Matrix<SCALAR, 3, 1>(
                      const Eigen::Matrix<double, 3, 1> &)>
                      u) {
  v_ = v;
  u_ = u;
}

private:
std::function<SCALAR(const Eigen::Matrix<double, 3, 1> &)> v_;
std::function<Eigen::Matrix<SCALAR, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
    u_;
SCALAR ana_sol;
};

}  // namespace debugging

}  // namespace projects::hldo_sphere

#endif  // THESIS_OPERATORS_WHITNEY_ONE_HODGE_LAPLACE_H
