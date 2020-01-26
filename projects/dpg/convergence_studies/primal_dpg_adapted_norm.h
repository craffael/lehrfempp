/**
 * @file
 * @headerfile projects/convergence_studies/primal_dpg_adapted_norm.h
 * @brief helper function to test convergence of the primal DPG method.
 * @author Philippe Peter
 * @date July 2019
 * @copyright MIT License
 */
#ifndef PROJECTS_DPG_CS_PRIMAL_DPG_AN_H
#define PROJECTS_DPG_CS_PRIMAL_DPG_AN_H

#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include "../dpg_element_matrix_provider.h"
#include "../dpg_element_vector_provider.h"
#include "../product_element_matrix_provider_builder.h"
#include "../product_element_vector_provider_builder.h"
#include "../product_fe_space.h"
#include "../product_fe_space_factory.h"

#include "../test/convection_diffusion_ell_bvp.h"

// utility type to return the reported energy norms.
// we return: number of dofs, H1-error, error estimator.
using energy_vector = std::vector<std::tuple<int, double, double>>;

namespace projects::dpg::test {
/**
 * @brief Examines the convergence behaviour of the primal DPG method for a
 * _pure_ Dirichlet convection-diffusion BVP, using and adpated test space norm.
 * @param reflevels the number of mesh refinements
 * @param solution exact solution
 * @param sol_grad  gradient of the exact solution
 * @param alpha scalar-valued diffusion coefficient of the BVP
 * @param beta advection field of the BVP
 * @param f source function of the BVP
 * @param g Dirichlet data of the BVP
 * @param degree_p polynomial degree of the u-part of the DPG approximation
 * @param enrichement value of \f$ \Delta p \f$, s.t. \f$ r  = p  + \Delta p \f$
 * is the polynomial degree of enriched test space
 * @param ref_el type of reference element for the tensor-product mesh
 * @return a vector containing for each refinement level the number of degrees
 * of freedom (f the trial space), the H1 error of the solution and the value of
 * the DPG error estimator.
 *
 * Tests the convergence behaviour of the primal DPG method on a tensor product
 * mesh based on uniform refinement. The method uses an adapted test space norm.
 */
template <typename SOLFUNC, typename SOLGRAD, typename ALPHAFUNC,
          typename BETAFUNC, typename FFUNC, typename GFUNC>
energy_vector
TestConververgencePrimalDPGAdaptedNormConvectionDiffusionDirichletBVP(
    size_type reflevels, SOLFUNC solution, SOLGRAD sol_grad, ALPHAFUNC alpha,
    BETAFUNC beta, FFUNC f, GFUNC g, int degree_p, int enrichement,
    lf::base::RefEl ref_el) {
  std::vector<std::tuple<int, double, double>> errnorms{};
  std::cout << "Approximation order: " << degree_p
            << ", Enrichement: " << enrichement << std::endl;

  // calculate degree for enriched spaces
  int degree_r = degree_p + enrichement;

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
    auto u = factory_trial.AddH1Component(degree_p);
    auto q_n = factory_trial.AddFluxComponent(degree_p - 1);
    auto fe_space_trial = factory_trial.Build();

    // construct test space:
    ProductUniformFESpaceFactory<double> factory_test(mesh_p);
    auto v = factory_test.AddL2Component(degree_r);
    auto fe_space_test = factory_test.Build();

    // wrap functions into a mesh function:
    auto alpha_mf = lf::mesh::utils::MeshFunctionGlobal(alpha);
    auto beta_mf = lf::mesh::utils::MeshFunctionGlobal(beta);
    auto f_mf = lf::mesh::utils::MeshFunctionGlobal(f);
    auto one_mf = lf::mesh::utils::MeshFunctionConstant(1.0);
    auto zero_mf =
        lf::mesh::utils::MeshFunctionConstant(Eigen::Vector2d(0.0, 0.0));

    // construct extended stiffness provider
    ProductElementMatrixProviderBuilder stiffness_builder(fe_space_trial,
                                                          fe_space_test);
    stiffness_builder.AddDiffusionElementMatrixProvider(u, v, alpha_mf);
    stiffness_builder.AddFluxElementMatrixProvider(q_n, v, -one_mf);
    stiffness_builder.AddConvectionElementMatrixProvider(u, v, zero_mf,
                                                         beta_mf);
    auto stiffness_provider = stiffness_builder.Build();

    // construct gram matrix provider
    auto alpha_squared_mf = lf::mesh::utils::MeshFunctionGlobal(
        [alpha](const Eigen::Vector2d& x) { return alpha(x) * alpha(x); });
    ProductElementMatrixProviderBuilder gramian_builder(fe_space_test,
                                                        fe_space_test);
    gramian_builder.AddDiffusionElementMatrixProvider(v, v, alpha_squared_mf);
    gramian_builder.AddReactionElementMatrixProvider(v, v, one_mf);
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
        u, q_n);

    // calculate full solution
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd sol_vec = solver.solve(phi);

    // extract u component:
    auto& dofh = fe_space_trial->LocGlobMap();
    auto fe_space_u = fe_space_trial->ComponentFESpace(u);
    Eigen::VectorXd sol_vec_u = sol_vec.head(dofh.NumDofs(u));

    // wrap finite element solution and exact solution into a mesh function:
    auto mf_fe_u = lf::uscalfe::MeshFunctionFE(fe_space_u, sol_vec_u);
    auto mf_fe_grad_u = lf::uscalfe::MeshFunctionGradFE(fe_space_u, sol_vec_u);
    auto mf_solution = lf::mesh::utils::MeshFunctionGlobal(solution);
    auto mf_solution_grad = lf::mesh::utils::MeshFunctionGlobal(sol_grad);

    // calculate the error in the H1 norm
    double H1_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(
        *mesh_p,
        lf::uscalfe::squaredNorm(mf_fe_u - mf_solution) +
            lf::uscalfe::squaredNorm(mf_fe_grad_u - mf_solution_grad),
        10));

    // evaluate error estimators:
    auto local_errors =
        ElementErrorEstimators(fe_space_trial->LocGlobMap(), stiffness_provider,
                               gramian_provider, rhs_provider, sol_vec);
    double posteriorError = EvalPosteriorError(mesh_p, local_errors);
    std::cout << "Level " << l << " ( " << dofh.NumDofs() << ")"
              << ", H1 error: " << H1_err << ", post:" << posteriorError
              << std::endl;
    errnorms.emplace_back(dofh.NumDofs(), H1_err, posteriorError);
  }
  return errnorms;
}
}  // namespace projects::dpg::test

#endif  // PROJECTS_DPG_CS_PRIMAL_DPG_AN_H
