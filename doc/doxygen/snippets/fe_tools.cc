/**
 * @file
 * @brief Doxygens snippets for functions in fe_tools.h
 * @author Raffael Casagrande
 * @date   2019-03-04 10:51:20
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe {

void integrateMeshFunction() {
  //! [integrateMeshFunction]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  // integrate the function sin(x)*cos(y) over the mesh using 5th-degree
  // quadrature rules
  auto mf = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  auto integral = IntegrateMeshFunction(*mesh, mf, 5);
  //! [integrateMeshFunction]
}

void integrateMeshFunction2() {
  //! [integrateMeshFunction2]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  // integrate the function sin(x)*cos(y) over the mesh using 5th-degree
  // quadrature rules
  auto mf = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return std::sin(x[0]) * std::cos(x[1]); });

  // select the quadrature rule explicitly for every element:
  auto integral = IntegrateMeshFunction(*mesh, mf, [](const mesh::Entity& e) {
    return quad::make_QuadRule(e.RefEl(), 5);
  });
  //! [integrateMeshFunction2]
}

void nodalProjection() {
  //! [nodalProjection]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  // project a first-degree polynomial onto a first order lagrange space and
  // make sure the representation is exact:
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);
  auto mf_linear = MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) { return 2 + 3 * x[0] + 4 * x[1]; });

  auto dof_vector = NodalProjection(*fe_space, mf_linear);
  auto mf_fe = MeshFunctionFE(fe_space, dof_vector);

  assert(IntegrateMeshFunction(*mesh, squaredNorm(mf_fe - mf_linear), 2) <
         1e-12);
  //! [nodalProjection]
}

void InitEssentialConditionFromFunction() {
  //! [InitEssentialConditionFromFunction]
  auto mesh_factory = std::make_unique<mesh::hybrid2d::MeshFactory>(2);
  auto gmsh_reader = io::GmshReader(std::move(mesh_factory), "mesh.msh");
  auto mesh = gmsh_reader.mesh();

  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh);

  // We want to solve the pde
  // - \laplace u = 0
  // u = cos(x)*sin(y) on the boundary

  // 1) Setup a mesh function that represents the prescribed values of u on the
  // boundary
  auto mf_boundary = MeshFunctionGlobal([&](const Eigen::Vector2d& x) {
    return std::cos(x[0]) * std::sin(x[1]);
  });

  // 2) determine the dofs on the boundary and to what value the should be set
  auto boundary_cond = InitEssentialConditionFromFunction(
      fe_space->LocGlobMap(),
      *fe_space->ShapeFunctionLayout(base::RefEl::kSegment()),
      base::PredicateTrue{}, mf_boundary);

  // 3) Assemble the stiffness matrix:
  assemble::COOMatrix<double> lhs(fe_space->LocGlobMap().NoDofs(),
                                  fe_space->LocGlobMap().NoDofs());
  auto mf_one = MeshFunctionConstant(1.);
  auto mf_zero = MeshFunctionConstant(0.);
  auto matrix_provider =
      ReactionDiffusionElementMatrixProvider(fe_space, mf_one, mf_zero);
  assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(),
                                  fe_space->LocGlobMap(), matrix_provider, lhs);

  // 4) Modify the system of equations such that the boundary values are
  // enforced
  Eigen::VectorXd rhs(fe_space->LocGlobMap().NoDofs());
  assemble::fix_flagged_solution_components(
      [&](unsigned int idx) { return boundary_cond[idx]; }, lhs, rhs);

  // 5) solve the problem:
  auto lhs_sparse = lhs.makeSparse();
  Eigen::SimplicialLDLT<decltype(lhs_sparse)> solver;
  solver.compute(lhs_sparse);
  auto x = solver.solve(rhs);

  //! [InitEssentialConditionFromFunction]
}

}  // namespace lf::uscalfe
