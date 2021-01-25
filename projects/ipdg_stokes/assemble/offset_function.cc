#include "offset_function.h"

#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/base/base.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/lagr_fe.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "utils.h"

namespace projects::ipdg_stokes::assemble {

PiecewiseBoundaryNormalJumpAssembler::PiecewiseBoundaryNormalJumpAssembler(
    std::shared_ptr<const lf::mesh::Mesh> mesh,
    const lf::mesh::utils::MeshDataSet<bool> &boundary)
    : mesh_(std::move(mesh)), boundary_(boundary) {}

Eigen::MatrixXd PiecewiseBoundaryNormalJumpAssembler::Eval(
    const lf::mesh::Entity &entity) const {
  // Use the vertex coordinates to compute the local normals on the edges
  const Eigen::Matrix<double, 2, 3> normals =
      projects::ipdg_stokes::mesh::computeOutwardNormals(entity);
  // Construct the basis functions from the curl of the standard hat functions
  const auto geom = entity.Geometry();
  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
  const Eigen::MatrixXd ref_grads =
      hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
          .transpose();
  const Eigen::MatrixXd J_inv_trans =
      geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));
  const Eigen::MatrixXd grads = J_inv_trans * ref_grads;
  Eigen::MatrixXd basis_funct(2, 3);
  basis_funct << grads.row(1), -grads.row(0);
  const auto nodes = entity.SubEntities(2);
  const auto edges = entity.SubEntities(1);
  // Count the number of boundary nodes and edges belonging to this triangle
  // to get the size of the element matrix
  int boundary_node_count = 0;
  int boundary_edge_count = 0;
  for (const auto node : nodes) {
    if (boundary_(*node)) {
      ++boundary_node_count;
    }
  }
  for (const auto edge : edges) {
    if (boundary_(*edge)) {
      ++boundary_edge_count;
    }
  }
  // Assemble the element matrix itself
  Eigen::MatrixXd elem_mat = Eigen::MatrixXd::Zero(
      std::max(1, boundary_edge_count), std::max(1, boundary_node_count));
  int col_idx = 0;
  for (int node_idx = 0; node_idx < 3; ++node_idx) {
    int row_idx = 0;
    const auto &node = *nodes[node_idx];
    if (boundary_(node)) {
      for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
        const auto &edge = *edges[edge_idx];
        if (boundary_(edge)) {
          elem_mat(row_idx, col_idx) =
              normals.col(edge_idx).transpose() * basis_funct.col(node_idx);
          ++row_idx;
        }
      }
      ++col_idx;
    }
  }
  return elem_mat;
}

Eigen::VectorXd createOffsetFunction(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::mesh::utils::MeshDataSet<bool> &boundary,
    const lf::assemble::DofHandler &dofh,
    const std::function<Eigen::Vector2d(const lf::mesh::Entity &)>
        &dirichlet_data,
    const Eigen::SparseMatrix<double> &A) {
  // Create a new dofhandler for the boundary nodes of the mesh
  lf::assemble::DynamicFEDofHandler boundary_node_dofh(
      mesh, [&](const lf::mesh::Entity &entity) {
        if (entity.RefEl() == lf::base::RefElType::kPoint && boundary(entity)) {
          return 1;
        }
        return 0;
      });
  // Create a new dofhandler for the boundary edges of the mesh
  lf::assemble::DynamicFEDofHandler boundary_edge_dofh(
      mesh, [&](const lf::mesh::Entity &entity) {
        if (entity.RefEl() == lf::base::RefElType::kSegment &&
            boundary(entity)) {
          return 1;
        }
        return 0;
      });
  // Assemble a matrix mapping the basis function coefficients on the boundary
  // nodes to the normal jumps over the boundary edges
  lf::assemble::COOMatrix<double> J(boundary_edge_dofh.NumDofs() + 1,
                                    boundary_node_dofh.NumDofs() + 1);
  PiecewiseBoundaryNormalJumpAssembler element_matrix_provider(mesh, boundary);
  lf::assemble::AssembleMatrixLocally(0, boundary_node_dofh, boundary_edge_dofh,
                                      element_matrix_provider, J);
  // Add the lagrange multiplier needed such that the average over all
  // coefficients is zero
  for (lf::base::size_type i = 0; i < boundary_node_dofh.NumDofs(); ++i) {
    J.AddToEntry(boundary_edge_dofh.NumDofs(), i, 1);
  }
  for (lf::base::size_type i = 0; i < boundary_edge_dofh.NumDofs(); ++i) {
    J.AddToEntry(i, boundary_node_dofh.NumDofs(), 1);
  }
  Eigen::SparseMatrix<double> Js = J.makeSparse();
  // Assemble the right hand side from the provided dirichlet data
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(boundary_edge_dofh.NumDofs() + 1);
  for (const auto cell : mesh->Entities(0)) {
    const Eigen::Matrix<double, 2, 3> normals =
        projects::ipdg_stokes::mesh::computeOutwardNormals(*cell);
    const auto edges = cell->SubEntities(1);
    for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
      const auto &edge = *edges[edge_idx];
      if (boundary(edge)) {
        rhs[boundary_edge_dofh.GlobalDofIndices(edge)[0]] =
            normals.col(edge_idx).transpose() * dirichlet_data(edge);
      }
    }
  }
  // Solve for the basis function coefficients of the offset function ordered
  // with the boundary dof handler
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(Js);
  const Eigen::VectorXd offset_function_local =
      solver.solve(rhs).head(boundary_node_dofh.NumDofs());
  // Compute the full offset function including the jump terms
  Eigen::VectorXd offset_function = Eigen::VectorXd::Zero(dofh.NumDofs());
  for (lf::assemble::size_type node_idx = 0;
       node_idx < boundary_node_dofh.NumDofs(); ++node_idx) {
    const auto &node = boundary_node_dofh.Entity(node_idx);
    offset_function[dofh.GlobalDofIndices(node)[0]] =
        offset_function_local[node_idx];
  }
  offset_function.tail(mesh->NumEntities(1)) =
      A.block(mesh->NumEntities(2), 0, mesh->NumEntities(1), A.cols()) *
      offset_function;
  return offset_function;
}

}  // end namespace projects::ipdg_stokes::assemble
