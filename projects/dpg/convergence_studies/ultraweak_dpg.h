/**
 * @file
 * @headerfile projects/dpg/convergence_studies/ultraweak_dpg.h
 * @brief helper function to test convergence of the ultraweak DGP method
 * @author Philippe Peter
 * @date July 2019
 * @copyright MIT License
 */
#ifndef PROJECTS_DPG_CS_ULTRAWEAK_DPG_H
#define PROJECTS_DPG_CS_ULTRAWEAK_DPG_H

#include <cmath>

#include <lf/base/base.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/fe/fe.h>

#include "../dpg_element_matrix_provider.h"
#include "../dpg_element_vector_provider.h"
#include "../dpg_tools.h"
#include "../product_element_matrix_provider_builder.h"
#include "../product_element_vector_provider_builder.h"
#include "../product_fe_space.h"
#include "../product_fe_space_factory.h"

#include "../test/convection_diffusion_ell_bvp.h"

// utility type to return the reported energy norms:
//[dofs, error in u-component,error in sigma-component, error estimate.]
using energy_vector_ultraweak =
    std::vector<std::tuple<int, double, double, double>>;

namespace projects::dpg::test {

/**
 * @brief Examines the convergence behaviour of the ultraweak DPG method for a
 * _pure_ Dirichlet convection-diffusion BVP
 * @param reflevels the number of mesh refinements
 * @param solution_u the \f$ u \f$-component of exact solution
 * @param solution_sigma \f$ \sigma \f$-component of the exact solution
 * @param alpha  scalar-valued diffusion coefficient of the BVP
 * @param beta  advection field of the BVP
 * @param f  source function of the BVP
 * @param g  Dirichlet Data of the BVP
 * @param degree polynoimial degree  of the \f$ u\f$-part of the DPG
 * approximation
 * @param enrichement value of \f$ \Delta p \f$, s.t. \f$ r  = p  + \Delta p \f$
 * is the polynomial degree of enriched test space
 * @param ref_el type of reference element for the tensor-product mesh
 * @return a vector, containing for each refinement level the number of degrees
 * of freedom (of the trial space), the \f$ L^2(\Omega) \f$ error in the u and
 * \$f\sigma\f$ component and the value of the DPG error estimator.
 *
 * Tests the convergence behaviour of the ultraweak DPG method on a tensor
 * product mesh based on uniform refinement.
 */
template <typename SOLFUNC_U, typename SOLFUNC_SIGMA, typename ALPHAFUNC,
          typename BETAFUNC, typename FFUNC, typename GFUNC>
energy_vector_ultraweak
TestConververgenceUltraWeakDPGConvectionDiffusionDirichletBVP(
    size_type reflevels, SOLFUNC_U solution_u, SOLFUNC_SIGMA solution_sigma,
    ALPHAFUNC alpha, BETAFUNC beta, FFUNC f, GFUNC g, int degree,
    int enrichement, lf::base::RefEl ref_el) {
  energy_vector_ultraweak errnorms{};

  // generate a tensor product mesh based on the provided reference element.
  std::shared_ptr<lf::mesh::Mesh> top_mesh;
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  switch (ref_el) {
    case lf::base::RefEl::kTria(): {
      // construct a triangular tensor product mesh
      lf::mesh::hybrid2d::TPTriagMeshBuilder builder(
          std::move(mesh_factory_ptr));
      builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
          .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
          .setNumXCells(2)
          .setNumYCells(2);
      top_mesh = builder.Build();
      break;
    }

    case lf::base::RefEl::kQuad(): {
      lf::mesh::hybrid2d::TPQuadMeshBuilder builder(
          std::move(mesh_factory_ptr));
      builder.setBottomLeftCorner(Eigen::Vector2d{0.0, 0.0})
          .setTopRightCorner(Eigen::Vector2d{1.0, 1.0})
          .setNumXCells(2)
          .setNumYCells(2);
      top_mesh = builder.Build();
      break;
    }
    default: {
      LF_ASSERT_MSG(false, "unsupported reference element:" << ref_el);
      break;
    }
  }

  // Generate a hierarchy of test meshes by uniform refinement
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(top_mesh,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  std::cout << multi_mesh;

  // get number of levels:
  size_type L = multi_mesh.NumLevels();

  // perform computations on all levels:
  for (size_type l = 0; l < L; ++l) {
    // retrive mesh for current level
    std::shared_ptr<const lf::mesh::Mesh> mesh_p{multi_mesh.getMesh(l)};

    // construct trial space:
    ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
    auto u = factory_trial.AddL2Component(degree);
    auto sigma_x = factory_trial.AddL2Component(degree);
    auto sigma_y = factory_trial.AddL2Component(degree);
    auto u_hat = factory_trial.AddTraceComponent(degree + 1);
    auto q_n = factory_trial.AddFluxComponent(degree);
    auto fe_space_trial = factory_trial.Build();

    // construct test space:
    ProductUniformFESpaceFactory<double> factory_test(mesh_p);
    auto v = factory_test.AddL2Component(degree + enrichement);
    auto tau_x = factory_test.AddL2Component(degree + enrichement);
    auto tau_y = factory_test.AddL2Component(degree + enrichement);
    auto fe_space_test = factory_test.Build();

    // wrap functions into a mesh function:
    auto alpha_mf = lf::mesh::utils::MeshFunctionGlobal(alpha);
    auto alpha_inf_mf = lf::mesh::utils::MeshFunctionGlobal(
        [&alpha](const Eigen::Vector2d& x) -> double {
          return 1.0 / alpha(x);
        });
    auto beta_mf = lf::mesh::utils::MeshFunctionGlobal(beta);
    auto f_mf = lf::mesh::utils::MeshFunctionGlobal(f);
    auto one_mf = lf::mesh::utils::MeshFunctionConstant(1.0);
    auto zero_mf =
        lf::mesh::utils::MeshFunctionConstant(Eigen::Vector2d(0.0, 0.0));
    auto x_selector_mf = lf::mesh::utils::MeshFunctionGlobal(
        [](const Eigen::Vector2d & /*x*/) -> Eigen::Matrix2d {
          return (Eigen::MatrixXd(2, 2) << 1.0, 0.0, 0.0, 0.0).finished();
        });
    auto y_selector_mf = lf::mesh::utils::MeshFunctionGlobal(
        [](const Eigen::Vector2d & /*x*/) -> Eigen::Matrix2d {
          return (Eigen::MatrixXd(2, 2) << 0.0, 0.0, 0.0, 1.0).finished();
        });

    auto x_selector_mf_vector =
        lf::mesh::utils::MeshFunctionConstant(Eigen::Vector2d(1.0, 0.0));
    auto y_selector_mf_vector =
        lf::mesh::utils::MeshFunctionConstant(Eigen::Vector2d(0.0, 1.0));

    // construct extended stiffness provider
    ProductElementMatrixProviderBuilder stiffness_builder(fe_space_trial,
                                                          fe_space_test);
    stiffness_builder.AddConvectionElementMatrixProvider(u, v, -beta_mf,
                                                         zero_mf);
    stiffness_builder.AddConvectionElementMatrixProvider(
        sigma_x, v, x_selector_mf_vector, zero_mf);
    stiffness_builder.AddConvectionElementMatrixProvider(
        sigma_y, v, y_selector_mf_vector, zero_mf);

    // stiffness_builder.AddTraceElementMatrixProvider(u_hat, v, beta_mf);
    stiffness_builder.AddFluxElementMatrixProvider(q_n, v, -one_mf);

    stiffness_builder.AddConvectionElementMatrixProvider(
        u, tau_x, x_selector_mf_vector, zero_mf);
    stiffness_builder.AddConvectionElementMatrixProvider(
        u, tau_y, y_selector_mf_vector, zero_mf);

    stiffness_builder.AddReactionElementMatrixProvider(sigma_x, tau_x,
                                                       alpha_inf_mf);
    stiffness_builder.AddReactionElementMatrixProvider(sigma_y, tau_y,
                                                       alpha_inf_mf);

    stiffness_builder.AddTraceElementMatrixProvider(u_hat, tau_x,
                                                    -x_selector_mf_vector);
    stiffness_builder.AddTraceElementMatrixProvider(u_hat, tau_y,
                                                    -y_selector_mf_vector);
    auto stiffness_provider = stiffness_builder.Build();

    // construct gram matrix provider
    ProductElementMatrixProviderBuilder gramian_builder(fe_space_test,
                                                        fe_space_test);
    gramian_builder.AddDiffusionElementMatrixProvider(v, v, one_mf);
    gramian_builder.AddReactionElementMatrixProvider(v, v, one_mf);

    gramian_builder.AddDiffusionElementMatrixProvider(tau_x, tau_x,
                                                      x_selector_mf);
    gramian_builder.AddDiffusionElementMatrixProvider(tau_y, tau_y,
                                                      y_selector_mf);

    gramian_builder.AddDiffusionElementMatrixProvider(
        tau_y, tau_x,
        lf::mesh::utils::MeshFunctionGlobal(
            [](const Eigen::Vector2d & /*x*/) -> Eigen::Matrix2d {
              return (Eigen::MatrixXd(2, 2) << 0.0, 1.0, 0.0, 0.0).finished();
            }));
    gramian_builder.AddDiffusionElementMatrixProvider(
        tau_x, tau_y,
        lf::mesh::utils::MeshFunctionGlobal(
            [](const Eigen::Vector2d & /*x*/) -> Eigen::Matrix2d {
              return (Eigen::MatrixXd(2, 2) << 0.0, 0.0, 1.0, 0.0).finished();
            }));

    gramian_builder.AddReactionElementMatrixProvider(tau_x, tau_x, one_mf);
    gramian_builder.AddReactionElementMatrixProvider(tau_y, tau_y, one_mf);
    auto gramian_provider = gramian_builder.Build();

    // construct the rhs provider:
    ProductElementVectorProviderBuilder rhs_builder(fe_space_test);
    rhs_builder.AddLoadElementVectorProvider(v, f_mf);
    auto rhs_provider = rhs_builder.Build();

    // initialize the dpg providers:
    auto element_matrix_provider =
        std::make_shared<DpgElementMatrixProvider<double>>(stiffness_provider,
                                                           gramian_provider);
    auto element_vector_provider =
        std::make_shared<DpgElementVectorProvider<double>>(
            rhs_provider, stiffness_provider, gramian_provider);

    // initialize the boundary value problem
    auto h = [](const Eigen::Vector2d & /*x*/) -> double { return 0.0; };
    auto dirichlet_selector = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p);
    auto bvp = std::make_shared<FullConvectionDiffusionBVP<
        decltype(alpha), decltype(beta), decltype(f), decltype(g), decltype(h),
        decltype(dirichlet_selector)>>(
        FullConvectionDiffusionBVP(alpha, beta, f, g, h,
                                   std::move(dirichlet_selector)));

    // Assemble full system:
    auto [A, phi] = ConvectionDiffusionDPGLinSys<double>(
        fe_space_trial, element_matrix_provider, element_vector_provider, bvp,
        u_hat, q_n);

    // calculate full solution
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd sol_vec = solver.solve(phi);

    // extract u component:
    auto& dofh = fe_space_trial->LocGlobMap();
    auto fe_space_u = fe_space_trial->ComponentFESpace(u);
    Eigen::VectorXd sol_vec_u =
        sol_vec.segment(dofh.Offset(u), dofh.NumDofs(u));

    // Extract sigma_x component:
    auto fe_space_sigma_x = fe_space_trial->ComponentFESpace(sigma_x);
    Eigen::VectorXd sol_vec_sigma_x =
        sol_vec.segment(dofh.Offset(sigma_x), dofh.NumDofs(sigma_x));

    // Extract sigma_y component:
    auto fe_space_sigma_y = fe_space_trial->ComponentFESpace(sigma_y);
    Eigen::VectorXd sol_vec_sigma_y =
        sol_vec.segment(dofh.Offset(sigma_y), dofh.NumDofs(sigma_y));

    // wrap finite element solution and exact solution into a mesh function:
    auto mf_fe_u = lf::fe::MeshFunctionFE(fe_space_u, sol_vec_u);
    auto mf_solution_u = lf::mesh::utils::MeshFunctionGlobal(solution_u);

    auto mf_fe_sigma_x =
        lf::fe::MeshFunctionFE(fe_space_sigma_x, sol_vec_sigma_x);
    auto mf_solution_sigma_x = lf::mesh::utils::MeshFunctionGlobal(
        [&solution_sigma](Eigen::Vector2d x) -> double {
          return solution_sigma(x)[0];
        });

    auto mf_fe_sigma_y =
        lf::fe::MeshFunctionFE(fe_space_sigma_y, sol_vec_sigma_y);
    auto mf_solution_sigma_y = lf::mesh::utils::MeshFunctionGlobal(
        [&solution_sigma](Eigen::Vector2d x) -> double {
          return solution_sigma(x)[1];
        });

    double L2err_u = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_fe_u - mf_solution_u), 10));

    double L2err_sigma_x = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_fe_sigma_x - mf_solution_sigma_x),
        10));
    double L2err_sigma_y = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p, lf::mesh::utils::squaredNorm(mf_fe_sigma_y - mf_solution_sigma_y),
        10));

    double L2err_sigma = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh_p,
        lf::mesh::utils::squaredNorm(mf_fe_sigma_x - mf_solution_sigma_x) +
            lf::mesh::utils::squaredNorm(mf_fe_sigma_y - mf_solution_sigma_y),
        10));

    auto local_errors =
        ElementErrorEstimators(fe_space_trial->LocGlobMap(), stiffness_provider,
                               gramian_provider, rhs_provider, sol_vec);
    double posteriorError = EvalPosteriorError(mesh_p, local_errors);
    std::cout << "Level " << l << "(" << dofh.NumDofs()
              << ") , L2 Error u: " << L2err_u
              << ", L2 error sigma: " << L2err_sigma << "(x:" << L2err_sigma_x
              << ", y: " << L2err_sigma_y << ")"
              << ", posterior error: " << posteriorError << std::endl;
    errnorms.emplace_back(dofh.NumDofs(), L2err_u, L2err_sigma, posteriorError);
  }
  return errnorms;
}

}  // namespace projects::dpg::test

#endif  // PROJECTS_DPG_CS_ULTRAWEAK_DPG_H
