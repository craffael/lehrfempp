/**
 * @file
 * @brief Solution of general second-order elliptic boundary value problem with
 * linear finite elements from a Gmsh generated mesh
 * @author Simon Meierhans
 * @date   January 2019
 * @copyright MIT License
 */

#include <fstream>
#include <iomanip>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/filesystem.hpp>

int main() {
  // abbreviations for types
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;

  // find path to mesh
  boost::filesystem::path here = __FILE__;
  auto mesh_path = here.parent_path() / "meshes/square.msh";

  // load mesh of square computational domain
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());

  // get pointer to mesh
  auto mesh = reader.mesh();

  // Count the number of edges with von Neumann boundary condition
  auto physical_entity_nr_neu = reader.PhysicalEntityName2Nr("neu");
  int num_neumann = 0;
  for (auto& e : mesh->Entities(1)) {
    if (reader.IsPhysicalEntity(e, physical_entity_nr_neu)) {
      ++num_neumann;
    }
  }
  std::cout << "boundary edges with von Neumann boundary condition: "
            << num_neumann << "\n";

  // Count the number of edges with Dirichlet boundary condition
  auto physical_entity_nr_dir = reader.PhysicalEntityName2Nr("dir");
  int num_dirichlet = 0;
  for (auto& e : mesh->Entities(1)) {
    if (reader.IsPhysicalEntity(e, physical_entity_nr_dir)) {
      ++num_dirichlet;
    }
  }
  std::cout << "boundary edges with Dirichlet boundary condition: "
            << num_dirichlet << "\n";

  // set up finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  // set up dof handler
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // Dimension of finite element space`
  const size_type N_dofs(dofh.NoDofs());

  // identity mesh function for very simple problem
  auto identity = [](Eigen::Vector2d x) -> double { return 1.; };
  lf::uscalfe::MeshFunctionGlobal mf_identity{identity};

  auto zero = [](Eigen::Vector2d x) -> double { return 0.; };
  lf::uscalfe::MeshFunctionGlobal mf_zero{zero};

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_identity), decltype(mf_identity)>
      elmat_builder(fe_space, mf_identity, mf_identity);

  // Invoke assembly on cells (co-dimension = 0)
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Right-hand side vector; has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // Initialize object taking care of local computations on all cells for the
  // source f. The source is the identity function
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_identity)>
      elvec_builder(fe_space, mf_identity);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  if (num_neumann > 0) {
    // Select von Neumann edges
    auto edge_sel_neu = [&reader, physical_entity_nr_neu](
                            const lf::mesh::Entity& edge) -> bool {
      return reader.IsPhysicalEntity(edge, physical_entity_nr_neu);
    };

    // Add contributions of von Neumann boundary conditions
    lf::uscalfe::ScalarLoadEdgeVectorProvider<double, decltype(mf_identity),
                                              decltype(edge_sel_neu)>
        elvec_builder_neu(fe_space, mf_identity, edge_sel_neu);
    AssembleVectorLocally(1, dofh, elvec_builder_neu, phi);
  }

  if (num_dirichlet > 0) {
    // Select Dirichlet edges
    auto edge_sel_dir = [&reader, physical_entity_nr_dir](
                            const lf::mesh::Entity& edge) -> bool {
      return reader.IsPhysicalEntity(edge, physical_entity_nr_dir);
    };

    // Obtain specification for shape functions on edges
    std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
        rsf_edge_p = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
    LF_ASSERT_MSG(rsf_edge_p != nullptr, "FE specification for edges missing");

    // Fetch flags and values for degrees of freedom located on Dirichlet
    // edges.
    // bd_flags strictly speaking would not be necessary here since only
    // boundary edges are flagged as 'dir' anyway. In other cases this might
    // however be necessary.
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
    auto ess_bdc_flags_values{lf::uscalfe::InitEssentialConditionFromFunction(
        dofh, *rsf_edge_p,
        [&edge_sel_dir, &bd_flags](const lf::mesh::Entity& edge) -> bool {
          return (bd_flags(edge) && edge_sel_dir(edge));
        },
        mf_zero)};
    // Eliminate Dirichlet dofs from linear system
    lf::assemble::fix_flagged_solution_components<double>(
        [&ess_bdc_flags_values](glb_idx_t gdof_idx) {
          return ess_bdc_flags_values[gdof_idx];
        },
        A, phi);
  }
  // Assembly completed: Convert COO matrix into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  // Compute Energy
  // Matrix in triplet format holding Stiffness matrix.
  if (N_dofs > 0) {
    lf::assemble::COOMatrix<double> Stiffness(N_dofs, N_dofs);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider<
        double, decltype(mf_identity), decltype(mf_zero)>
        stiffness_mat_builder(fe_space, mf_identity, mf_zero);
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stiffness_mat_builder,
                                        Stiffness);
    Eigen::SparseMatrix<double> Stiffness_crs = Stiffness.makeSparse();

    // Matrix in triplet format holding Mass matrix.
    lf::assemble::COOMatrix<double> Mass(N_dofs, N_dofs);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider<
        double, decltype(mf_zero), decltype(mf_identity)>
        mass_mat_builder(fe_space, mf_zero, mf_identity);
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, mass_mat_builder, Mass);
    Eigen::SparseMatrix<double> Mass_crs = Mass.makeSparse();

    // energy_stiffness_sq = 1' A 1
    double energy_stiffness_sq =
        Eigen::VectorXd::Constant(N_dofs, 1.0).transpose() *
        (Stiffness_crs * Eigen::VectorXd::Constant(N_dofs, 1.0));
    // energy_mass_sq = \mu' M \mu
    double energy_mass_sq = sol_vec.transpose() * Mass_crs * sol_vec;

    std::cout << "Computed Energy Norm: "
              << std::sqrt(energy_stiffness_sq + energy_mass_sq) << "\n";
  }
}
