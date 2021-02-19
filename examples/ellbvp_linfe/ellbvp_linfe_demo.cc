/**
 * @file
 * @brief Solution of general second-order elliptic boundary value problem with
 * linear finite elements
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <fstream>
#include <iomanip>

int main(int /*argc*/, const char** /*argv*/) {
  std::cout << "\t LehrFEM++ Demonstration Code " << std::endl;
  std::cout << "\t Solution of general second-order elliptic\n"
            << "\t boundary value problem by means of linear\n"
            << "\t Lagrangian finite element discretization" << std::endl;

  // abbreviations for types
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;

  // ======================================================================
  // Prep stage: provide all coefficient functions mainly through lambda
  //             functions and derived MeshFunctions.
  // ======================================================================

  // Coefficients:

  // 2x2 diffusion tensor A(x)
  auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix<double, 2, 2> {
    return (Eigen::Matrix<double, 2, 2>() << (3.0 + x[1]), x[0], x[0],
            (2.0 + x[0]))
        .finished();
  };
  // Wrap diffusion coefficient into a MeshFunction
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};

  // Scalar valued reaction coefficient c
  auto gamma = [](Eigen::Vector2d x) -> double {
    return (x[0] * x[0] + x[1] * x[1]);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_gamma{gamma};

  // Scalar valued impedance coefficient
  auto eta = [](Eigen::Vector2d x) -> double { return (1.0 + x[0] + x[1]); };
  lf::mesh::utils::MeshFunctionGlobal mf_eta{eta};

  /* SAM_LISTING_BEGIN_1 */
  // Exact solution u
  auto u = [](Eigen::Vector2d x) -> double {
    return std::log(x[0] * x[0] + x[1] + 1.0);
  };
  // Has to be wrapped into a mesh function for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_u{u};

  // Gradient of exact solution
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    double den = x[0] * x[0] + x[1] + 1.0;
    return ((Eigen::Vector2d() << 2.0 * x[0], 1.0).finished()) / den;
  };
  // Convert into mesh function to use for error computation
  lf::mesh::utils::MeshFunctionGlobal mf_grad_u{grad_u};
  /* SAM_LISTING_END_1 */

  // Right-hand side source function f
  auto f = [&gamma, &u](Eigen::Vector2d x) -> double {
    const double den = x[0] * x[0] + x[1] + 1.0;
    const double num =
        -x[0] * x[0] * (2 * x[1] + 9) - x[0] + 2 * x[1] * x[1] + 9 * x[1] + 5;
    return (-num / (den * den) + gamma(x) * u(x));
  };
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // Boundary conditions and predicates for different boundary parts
  // The predicates are named edge_sel_* and when invoked for an entity of edge
  // type return a boolean value
  // Auxiliary vector: Normal vector for Neumann boundary
  const Eigen::Vector2d n_neu =
      ((Eigen::Vector2d() << 1.0, 0.8).finished()).normalized();
  // Auxiliary vector; Normal vector for impedance boundary
  const Eigen::Vector2d n_imp =
      ((Eigen::Vector2d() << -1.0, 0.2).finished()).normalized();
  // Dirichlet data borrowed from the known exact solution.
  auto g = [&u](const Eigen::Vector2d& x) -> double { return u(x); };
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};

  // Predicates and selectors for boundary conditions
  // Predicates for selecting edges on the Neumann boundary
  auto neu_sel = [](Eigen::Vector2d x) -> bool {
    return (x[1] / 1.25 + x[0] - 1 > -1.0E-7);
  };
  lf::refinement::EntityCenterPositionSelector<
      std::function<bool(const Eigen::Vector2d&)>>
      edge_sel_neu{neu_sel};
  // Predicates for selecting edges on the Dirichlet boundary
  std::function<bool(const Eigen::Vector2d&)> dir_sel =
      [](const Eigen::Vector2d& x) -> bool { return (x[1] < 1.0E-7); };
  lf::refinement::EntityCenterPositionSelector<
      std::function<bool(const Eigen::Vector2d&)>>
      edge_sel_dir{dir_sel};
  // Predicates for selecting edges on the impedance boundary
  auto imp_sel = [](const Eigen::Vector2d& x) -> bool {
    return (x[1] - 5.0 * x[0] > -1.0E-7);
  };
  // Set up a predicate: this object is a functor taking an entity and returning
  // 'true' or 'false' depending on the location of the center of gravity of
  // that entity. In this demo code this is the way how to tell, which type of
  // boundary conditions prevails at different parts of the boundary.
  lf::refinement::EntityCenterPositionSelector edge_sel_imp{imp_sel};

  // Source term h on the right-hand side of Neumann and impedance boundary
  // conditions
  auto h = [&](const Eigen::Vector2d& x) -> double {
    if (imp_sel(x)) {
      // Impedance boundary
      return ((alpha(x) * grad_u(x)).dot(n_imp) + eta(x) * u(x));
    }
    if (neu_sel(x)) {
      // Neumann boundary
      return ((alpha(x) * grad_u(x)).dot(n_neu));
    }
    LF_ASSERT_MSG(false, "h called for Dirichlet edge!");
    return 0.0;
  };
  lf::mesh::utils::MeshFunctionGlobal mf_h{h};

  // ======================================================================
  // Stage I: Definition of computational domain through coarsest mesh
  // Since this example relies on a manufactured solution tied to a particular
  // domain, using a hard-wired mesh is justified. Another example will address
  // solving a boundary value problem on a mesh read from file.
  // ======================================================================

  // The following code also illustrates the role of a MeshFactory
  // Create helper object: mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Generate nodes of the mesh
  // clang-format off
  std::array<std::array<double, 2>, 8> node_coord{
  std::array<double, 2>({0   , 0    }),
  std::array<double, 2>({1   , 0    }),
  std::array<double, 2>({0.2 , 1.0  }),
  std::array<double, 2>({0.5 , 0.0  }),
  std::array<double, 2>({0.6 , 0.5  }),
  std::array<double, 2>({0.1 , 0.5  }),
  std::array<double, 2>({0.15, 0.75 }),
  std::array<double, 2>({0.4 , 0.75 })};
  // clang-format on
  // Add nodes to the mesh via the MeshFactory object
  for (const auto& node : node_coord) {
    mesh_factory_ptr->AddPoint(coord_t({node[0], node[1]}));
  }
  // Add plain triangles to the mesh, defined by their vertex nodes.
  // Since no particular geometry is specified, the triangles are assumed to
  // have straght edges.
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({3, 1, 4}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({7, 6, 2}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Create a general quadrilateral with straight edges.
  std::array<size_type, 4> quad_nodes{5, 4, 7, 6};
  Eigen::Matrix<double, 2, 4> quad_coord;
  for (int n_pt = 0; n_pt < 4; ++n_pt) {
    quad_coord(0, n_pt) = node_coord[quad_nodes[n_pt]][0];
    quad_coord(1, n_pt) = node_coord[quad_nodes[n_pt]][1];
  }
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(), std::vector<size_type>({5, 4, 7, 6}),
      std::make_unique<lf::geometry::QuadO1>(quad_coord));
  // Create parallelogram.
  std::array<size_type, 4> parg_nodes{0, 3, 4, 5};
  for (int n_pt = 0; n_pt < 4; ++n_pt) {
    quad_coord(0, n_pt) = node_coord[parg_nodes[n_pt]][0];
    quad_coord(1, n_pt) = node_coord[parg_nodes[n_pt]][1];
  }
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(), std::vector<size_type>({0, 3, 4, 5}),
      std::make_unique<lf::geometry::Parallelogram>(quad_coord));

  // Get a pointer to the coarsest mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  // Print information about the coarsest mesh
  std::cout << "\t Coarsest mesh for demonstration run\n";
  lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

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
  // Number of levels
  size_type L = multi_mesh.NumLevels();

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

    // Preprocessing: count number of edges with different boundary conditions
    size_type no_Dirichlet_edges = 0;
    size_type no_Neumann_edges = 0;
    size_type no_impedance_edges = 0;
    // Obtain an array of boolean flags for the edges of the mesh, 'true'
    // indicates that the edge lies on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
    // Traverse the edges of the mesh and check their boundary flags and the
    // type of boundary condition
    for (const lf::mesh::Entity* edge : mesh.Entities(1)) {
      if (bd_flags(*edge)) {
        if (edge_sel_dir(*edge)) {
          no_Dirichlet_edges++;
        } else if (edge_sel_imp(*edge)) {
          no_impedance_edges++;
        } else {
          no_Neumann_edges++;
        }
      }
    }
    // Dimension of finite element space`
    const size_type N_dofs(dofh.NumDofs());

    // Verbose output
    std::cout << "Computations on level " << level
              << " (#Dir_ed =  " << no_Dirichlet_edges
              << ", #Neu_ed = " << no_Neumann_edges
              << ", #imp_ed = " << no_impedance_edges << "), #dof = " << N_dofs
              << std::endl;

    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

    // ----------------------------------------------------------------------
    // III: Assemble finite element Galerkin matrix
    // First the volume part for the bilinear form
    // Initialize object taking care of local computations. No selection of a
    // subset of cells is specified in this demonstration: assembly will cover
    // all cells.
    lf::uscalfe::ReactionDiffusionElementMatrixProvider<
        double, decltype(mf_alpha), decltype(mf_gamma)>
        elmat_builder(fe_space, mf_alpha, mf_gamma);
    // Invoke assembly on cells (co-dimension = 0 as first argument)
    // Information about the mesh and the local-to-global map is passed through
    // a Dofhandler object, argument 'dofh'. This function call adds triplets to
    // the internal COO-format representation of the sparse matrix A.
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
    // Update with potential contributions from edges (boundary part of bilinear
    // form due to impedance boundary conditions!). To that end initialize an
    // object taking care of local computations on edges. A predicate, that is,
    // a functor returning boolean values, must ensure that computations are
    // confined to edges on the impedance boundary! This predicate is defined
    // next using the location-based predicate for impedance edges.
    auto imp_edge_sel = [&bd_flags,
                         &edge_sel_imp](const lf::mesh::Entity& edge) -> bool {
      return (bd_flags(edge) && edge_sel_imp(edge));
    };
    lf::uscalfe::MassEdgeMatrixProvider<double, decltype(mf_eta),
                                        decltype(imp_edge_sel)>
        edgemat_builder(fe_space, mf_eta, imp_edge_sel);
    // Invoke assembly on edges by specifying co-dimension = 1
    lf::assemble::AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, A);

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

    // Add contributions from Neumann and impedance edges. Note that that
    // right-hand side contribution of these two boundary conditions are given
    // by the same expression involving the source function h.
    if ((no_Neumann_edges > 0) || (no_impedance_edges > 0)) {
      auto edge_sel = [&bd_flags, &edge_sel_neu,
                       &edge_sel_imp](const lf::mesh::Entity& edge) -> bool {
        return (bd_flags(edge) && (edge_sel_neu(edge) || edge_sel_imp(edge)));
      };
      // Object taking care of local computations. A predicate selects the edges
      // to be processed
      lf::uscalfe::ScalarLoadEdgeVectorProvider<double, decltype(mf_h),
                                                decltype(edge_sel)>
          elvec_builder_neu(fe_space, mf_h, edge_sel);
      // Invoke assembly on edges (codim == 1), update vector
      AssembleVectorLocally(1, dofh, elvec_builder_neu, phi);
    }

    // ----------------------------------------------------------------------
    // III: Fixing solution components according to essential (Dirichlet)
    // boundary conditions
    if (no_Dirichlet_edges > 0) {
      // Obtain specification for shape functions on edges
      const auto* rsf_edge_p =
          fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
      LF_ASSERT_MSG(rsf_edge_p != nullptr,
                    "FE specification for edges missing");

      // Fetch flags and values for degrees of freedom located on Dirichlet
      // edges.
      auto ess_bdc_flags_values{lf::fe::InitEssentialConditionFromFunction(
          *fe_space,
          [&edge_sel_dir, &bd_flags](const lf::mesh::Entity& edge) -> bool {
            return (bd_flags(edge) && edge_sel_dir(edge));
          },
          mf_g)};
      // Eliminate Dirichlet dofs from linear system
      lf::assemble::FixFlaggedSolutionComponents<double>(
          [&ess_bdc_flags_values](glb_idx_t gdof_idx) {
            return ess_bdc_flags_values[gdof_idx];
          },
          A, phi);
    }

    /* SAM_LISTING_BEGIN_2 */
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
    /* SAM_LISTING_END_2 */
    errs.emplace_back(N_dofs, L2err, H1serr);
  }

  // Output table of errors to file and terminal
  std::ofstream out_file("errors.txt");
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
