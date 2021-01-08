#include "build_system_matrix.h"

#include <vector>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/lagr_fe.h>

#include "offset_function.h"
#include "piecewise_const_element_matrix_provider.h"
#include "piecewise_const_element_vector_provider.h"
#include "solution_to_mesh_data_set.h"
#include "utils.h"

namespace projects::ipdg_stokes::assemble {

std::tuple<lf::assemble::COOMatrix<double>, Eigen::VectorXd>
buildSystemMatrixNoFlow(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh,
    const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> &f,
    const std::function<Eigen::Vector2d(const lf::mesh::Entity &)>
        &dirichlet_data,
    double sigma, const lf::quad::QuadRule &quadrule, bool modified_penalty) {
  // Compute the boundary elements
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  // Compute the dirichlet data
  auto dirichlet = lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d>(mesh, 1);
  for (const auto *ep : mesh->Entities(1)) {
    dirichlet(*ep) = dirichlet_data(*ep);
  }
  // Assemble the system matrix
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  const auto elem_mat_builder =
      projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider(
          sigma, boundary, modified_penalty);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_builder, A);
  // Assemble the right hand side
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  const auto elem_vec_builder =
      projects::ipdg_stokes::assemble::PiecewiseConstElementVectorProvider(
          sigma, f, quadrule, boundary, dirichlet);
  lf::assemble::AssembleVectorLocally(0, dofh, elem_vec_builder, rhs);

  // Enforce no-flow boundary conditions
  auto selector = [&](lf::base::size_type idx) -> std::pair<bool, double> {
    const auto &e = dofh.Entity(idx);
    return {e.RefEl() == lf::base::RefElType::kPoint && boundary(e), 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A, rhs);

  // Return the computed LSE
  return {A, rhs};
}

std::tuple<lf::assemble::COOMatrix<double>, Eigen::VectorXd, Eigen::VectorXd>
buildSystemMatrixInOutFlow(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh,
    const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> &f,
    const std::function<Eigen::Vector2d(const lf::mesh::Entity &)>
        &dirichlet_data,
    double sigma, const lf::quad::QuadRule &quadrule, bool modified_penalty) {
  // Compute the boundary elements
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  // Compute the dirichlet data
  auto dirichlet = lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d>(mesh, 1);
  for (const auto *ep : mesh->Entities(1)) {
    dirichlet(*ep) = dirichlet_data(*ep);
  }
  // Assemble the system matrix
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  const auto elem_mat_builder =
      projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider(
          sigma, boundary, modified_penalty);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_builder, A);
  // Assemble the right hand side
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  const auto elem_vec_builder =
      projects::ipdg_stokes::assemble::PiecewiseConstElementVectorProvider(
          sigma, f, quadrule, boundary, dirichlet);
  lf::assemble::AssembleVectorLocally(0, dofh, elem_vec_builder, rhs);

  // Compute the offset function
  const Eigen::VectorXd offset_function = createOffsetFunction(
      mesh, boundary, dofh, dirichlet_data, A.makeSparse());

  // Apply offset function technique to the LSE
  rhs -= A.makeSparse().block(0, 0, dofh.NumDofs(), mesh->NumEntities(2)) *
         offset_function.head(mesh->NumEntities(2));
  auto selector = [&](lf::base::size_type idx) -> std::pair<bool, double> {
    const auto &e = dofh.Entity(idx);
    return {e.RefEl() == lf::base::RefElType::kPoint && boundary(e), 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A, rhs);

  // Return the computed LSE
  return {A, rhs, offset_function};
}

}  // end namespace projects::ipdg_stokes::assemble
