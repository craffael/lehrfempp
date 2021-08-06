/**
 * @file
 * @brief Solution of general second-order elliptic boundary value problem with
 * linear finite elements
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   Modified version for simple Dirichlet problem
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <fstream>
#include <iomanip>

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "\t LehrFEM++ Demonstration Code" << std::endl;
  std::cout << "\t Solution of general second-order elliptic\n"
            << "\t homogeneous Dirichlet problem by means of linear\n"
            << "\t Lagrangian finite element discretization" << std::endl;

  // abbreviations for types
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;

  // ======================================================================
  // Prep stage: provide all coefficient functions mainly through lambda
  //             functions and derived MeshFunctions.
  // ======================================================================

  // Coefficients:
 auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
   return Eigen::Matrix<double, 2, 2>::Identity();
 };
 auto u = [](Eigen::Vector2d x) -> double {
   return (std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]));
 };
 auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
   return M_PI *
          ((Eigen::Vector2d() << std::cos(M_PI * x(0)) * std::sin(M_PI * x(1)),
            std::sin(M_PI * x(0)) * std::cos(M_PI * x(1)))
               .finished());
 };
 auto f = [&u](Eigen::Vector2d x) -> double { return (2.0*M_PI * M_PI * u(x)); };
 /*
  auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    return Eigen::Matrix<double, 2, 2>::Identity();
  };
  auto u = [](Eigen::Vector2d x) -> double {
    return (x[0] * (1.0 - x[0]) * x[1] * (1.0 - x[1]));
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return ((Eigen::Vector2d() << (1 - 2 * x[0]) * x[1] * (1 - x[1]),
             x[0] * (1 - x[0]) * (1 - 2 * x[1]))
                .finished());
  };
  auto f = [](Eigen::Vector2d x) -> double {
    return 2 * (x[1] * (1 - x[1]) + x[0] * (1 - x[0]));
    }; */

  // Wrap diffusion coefficient into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};
  // Has to be wrapped into a mesh function for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  // Convert into mesh function to use for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // ======================================================================
  // Stage I: Definition of computational domain through coarsest mesh
  // Since this example relies on a manufactured solution tied to a particular
  // domain, using a hard-wired mesh is justified. Another example will address
  // solving a boundary value problem on a mesh read from file.
  // ======================================================================

  // A triangular mesh of the unit square
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);

  // ======================================================================
  // Stage II: Ask LehrFEM++ to create a hierarchy of nested meshes
  // ======================================================================

  // Obtain a pointer to a hierarchy of nested meshes
  const int reflevels = 7;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Ouput information about hierarchy of nested meshes
  std::cout << "\t Sequence of nested meshes used in demo code\n";
  multi_mesh.PrintInfo(std::cout);
  size_type L = multi_mesh.NumLevels();  // Number of levels

  // Vector for keeping error norms
  std::vector<std::tuple<size_type, double, double>> errs{};
  // LEVEL LOOP: Do computations on all levels
  for (size_type level = 0; level < L; ++level) {
    mesh_p = multi_mesh.getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Reference to current mesh
    const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

    // Dimension of finite element space`
    const size_type N_dofs(dofh.NumDofs());

    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

    // ----------------------------------------------------------------------
    // III: Assemble finite element Galerkin matrix
    // Helper object for computation of element matrices
    lf::fe::DiffusionElementMatrixProvider<double, decltype(mf_alpha)>
        elmat_builder(fe_space, mf_alpha);
    // Invoke assembly on cells (co-dimension = 0 as first argument)
    // Information about the mesh and the local-to-global map is passed through
    // a Dofhandler object, argument 'dofh'. This function call adds triplets to
    // the internal COO-format representation of the sparse matrix A.
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

    // ----------------------------------------------------------------------
    // IV: Right-hand side vector; has to be set to zero initially
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
    phi.setZero();
    // Assemble volume part of right-hand side vector depending on the source
    // function f.
    // Initialize object taking care of local computations on all cells.
    lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
        elvec_builder(fe_space, mf_f);
    // Invoke assembly on cells (codim == 0)
    AssembleVectorLocally(0, dofh, elvec_builder, phi);

    // ----------------------------------------------------------------------
    // III: Fixing solution components according to essential (Dirichlet)
    // boundary conditions
    // Obtain an array of boolean flags for the edges of the mesh, 'true'
    // indicates that the edge lies on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 2)};
    // Elimination of degrees of freedom on the boundary
    lf::assemble::FixFlaggedSolutionComponents<double>(
        [&bd_flags, &dofh](glb_idx_t gdof_idx) -> std::pair<bool, double> {
          const lf::mesh::Entity& node{dofh.Entity(gdof_idx)};
          return (bd_flags(node) ? std::make_pair(true, 0.0)
                                 : std::make_pair(false, 0.0));
        },
        A, phi);

    // Assembly completed: Convert COO matrix A into CRS format using Eigen's
    // internal conversion routines.
    Eigen::SparseMatrix<double> A_crs = A.makeSparse();

    // Solve linear system using Eigen's sparse direct elimination
    // Examine return status of solver in case the matrix is singular
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_crs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
    Eigen::VectorXd sol_vec = solver.solve(phi);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

    // Postprocessing: Compute error norms
    // create mesh functions representing solution / gradient of solution
    const lf::fe::MeshFunctionFE mf_sol(fe_space, sol_vec);
    const lf::fe::MeshFunctionGradFE mf_grad_sol(fe_space, sol_vec);
    // compute errors with 3rd order quadrature rules, which is sufficient for
    // piecewise linear finite elements
    double L2err =  // NOLINT
        std::sqrt(lf::fe::IntegrateMeshFunction(
            mesh, lf::mesh::utils::squaredNorm(mf_sol - mf_u), 2));
    double H1serr = std::sqrt(lf::fe::IntegrateMeshFunction(  // NOLINT
        mesh, lf::mesh::utils::squaredNorm(mf_grad_sol - mf_grad_u), 2));
    errs.emplace_back(N_dofs, L2err, H1serr);
  }

  // Output table of errors to file and terminal
  std::ofstream out_file("errors.txt");
  std::cout << "\t Table of error norms" << std::endl;
  std::cout << std::left << std::setw(10) << "N" << std::right << std::setw(16)
            << "L2 error" << std::setw(16) << "H1 error" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for (const auto& err : errs) {
    auto [N, l2err, h1serr] = err;
    out_file << std::left << std::setw(10) << N << std::left << std::setw(16)
             << l2err << std::setw(16) << h1serr << std::endl;
    std::cout << std::left << std::setw(10) << N << std::left << std::setw(16)
              << l2err << std::setw(16) << h1serr << std::endl;
  }

  return 0;
}
