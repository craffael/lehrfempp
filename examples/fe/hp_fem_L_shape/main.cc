/**
 * @file
 * @brief Example that shows how a singularity on an L-shaped domain can be
 * approximated with successive hp-fem refinements.
 * @author Raffael Casagrande
 * @date   2021-01-20 03:57:44
 * @copyright MIT License
 */

#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>

#include <filesystem>

int main() {
  // 1) Load the mesh:
  std::filesystem::path path(__FILE__);
  path = path.parent_path() / "L_shape.msh";
  lf::io::GmshReader reader(
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2), path.string());

  // 2) define the exact, analytic solution:
  auto exact_solution =
      lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
        double r = x.norm();
        double theta = std::atan2(x.y(), x.x());
        return std::pow(r, 2. / 3.) *
               std::sin(2. / 3. * (theta + lf::base::kPi / 2.));
      });
  auto diffusion_coefficient = lf::mesh::utils::MeshFunctionConstant(1.);

  // 3) Create mesh hierarchy and loop over levels:
  lf::refinement::MeshHierarchy mesh_hierarchy(
      reader.mesh(), std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2));
  const lf::base::size_type max_level = 9;
  Eigen::VectorXd errors(max_level + 1);
  Eigen::VectorXi num_dofs(max_level + 1);
  for (lf::base::size_type level = 0; level <= max_level; ++level) {
    auto mesh = mesh_hierarchy.getMesh(level);

    // create a MeshDataSet that stores for every point the adjacent cells:
    lf::mesh::utils::CodimMeshDataSet<std::vector<const lf::mesh::Entity*>>
        adjacent_cells(mesh, 2);
    for (const auto& ep : mesh->Entities(0)) {
      for (const auto& sp : ep->SubEntities(2)) {
        adjacent_cells(*sp).push_back(ep);
      }
    }

    // find the point entity at the origin (0,0):
    Eigen::Matrix<double, 0, 1> zero;
    const auto* const origin = std::find_if(
        mesh->Entities(2).begin(), mesh->Entities(2).end(), [&](auto ep) {
          auto coord = ep->Geometry()->Global(zero);
          return coord.norm() < 1e-10;
        });

    // 4) create a hp-approximation space where the degree of the elements
    // adjacent to the origin have degree 1. From there the degree is then
    // linearly increasing to level+1.
    lf::mesh::utils::AllCodimMeshDataSet<unsigned> degrees(mesh, level + 1);
    // -> mark the degree of all elements around the origin to be l1:
    for (auto& ep : adjacent_cells(**origin)) {
      degrees(*ep) = 1;
    }
    for (unsigned i = 0; i < level; ++i) {
      for (const auto& ep : mesh->Entities(0)) {
        if (degrees(*ep) == 1 + i) {
          // mark all adjacent cells to have one level more:
          for (const auto& sp : ep->SubEntities(2)) {
            for (const auto& neighbor : adjacent_cells(*sp)) {
              degrees(*neighbor) = std::min(degrees(*neighbor), 2 + i);
            }
          }
        }
      }
    }
    // set the edges to the minimum degree of the two adjacent elements:
    for (const auto* edge_p : mesh->Entities(1)) {
      degrees(*edge_p) = 1000;
    }
    for (const auto* ep : mesh->Entities(0)) {
      auto element_degree = degrees(*ep);
      for (const auto* edge_p : ep->SubEntities(1)) {
        degrees(*edge_p) = std::min(degrees(*edge_p), element_degree);
      }
    }
    auto fes = std::make_shared<lf::fe::HierarchicScalarFESpace<double>>(
        mesh, degrees);

    // 5) solve laplace problem with inhomogeneous dirichlet boundary conditions
    lf::base::size_type num_dof = fes->LocGlobMap().NumDofs();
    auto on_boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1);
    lf::fe::DiffusionElementMatrixProvider diffusion_emp(fes,
                                                         diffusion_coefficient);
    lf::assemble::COOMatrix<double> lhs_coo(num_dof, num_dof);
    AssembleMatrixLocally(0, fes->LocGlobMap(), fes->LocGlobMap(),
                          diffusion_emp, lhs_coo);

    Eigen::VectorXd rhs(num_dof);
    rhs.setZero();

    auto bnd_cond = lf::fe::InitEssentialConditionFromFunction(
        *fes, [&](const lf::mesh::Entity& e) { return on_boundary(e); },
        exact_solution);
    lf::assemble::FixFlaggedSolutionComponents(
        [&](auto i) { return bnd_cond[i]; }, lhs_coo, rhs);
    auto lhs = lhs_coo.makeSparse();

    std::cout << "Solving level " << level << " with " << num_dof << " dofs";
    num_dofs(level) = num_dof;
    Eigen::SimplicialLLT<decltype(lhs)> solver(lhs);
    LF_VERIFY_MSG(solver.info() == Eigen::Success,
                  "Sparse Cholesky decomposition was not successful.");
    auto solution_vec = solver.solve(rhs).eval();
    LF_VERIFY_MSG(solver.info() == Eigen::Success,
                  "problem while applying sparse cholesky decomposition.");

    // 6) Compute L2 error to analytical solution
    auto solution_mf = lf::fe::MeshFunctionFE(fes, solution_vec);

    auto L2error = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh, squaredNorm(solution_mf - exact_solution),
        static_cast<int>(2 * (level + 5))));
    errors(level) = L2error;
    std::cout << ", L2 error: " << L2error << std::endl;

    // lf::fe::DiffusionElementMatrixProvider laplacian()

    lf::io::VtkWriter vtk(mesh, "L_shape_" + std::to_string(level) + ".vtk", 0,
                          level + 1);
    vtk.WriteCellData("Degree", degrees);
    vtk.WritePointData("solution", solution_mf);

    if (level < max_level) {
      // refine the mesh towards the singularity:
      lf::mesh::utils::CodimMeshDataSet<bool> marked_edges(mesh, 1, false);
      for (const auto* ep : mesh->Entities(0)) {
        if (degrees(*ep) == 1) {
          // mark all adjacent edges:
          for (const auto* edge_p : ep->SubEntities(1)) {
            marked_edges(*edge_p) = true;
          }
        }
      }
      mesh_hierarchy.MarkEdges(
          [&](const lf::mesh::Mesh& /*unused*/, const lf::mesh::Entity& edge) {
            return marked_edges(edge);
          });
      mesh_hierarchy.RefineMarked();
    }
  }

  fmt::print("# dofs      : {}\n", num_dofs.transpose());
  fmt::print("max p degree: {}\n",
             Eigen::RowVectorXi::LinSpaced(max_level + 1, 1, max_level + 1));
  fmt::print("L2 errors   : {}\n", errors.transpose());
}
