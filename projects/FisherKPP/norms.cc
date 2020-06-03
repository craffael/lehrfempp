/** @file norms.cc
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 05.05.20
 *  @copyright ETH Zurich
 */

#include "norms.h"

#include <cmath>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

namespace FisherKPP {


double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p) {
  
  double mesh_size = 0.0;

  // Find maximal edge length
  double edge_length;
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    // Compute the length of the edge
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }

  return mesh_size;
} 

/*
 * Converts a large vector on  a grid to a smaller vector correponding to a sub-grid.
 * mu vector of function values on a spacial grid of size N_large
 * N divides the size N_large of mu
 * returns a vector mu_sub of size N, that represents mu on a sub-grid of size N
 */
Eigen::VectorXd reduce(const Eigen::VectorXd &mu, unsigned int N) {
  Eigen::VectorXd mu_sub(N);
  int fraction = mu.size() / N;
  for (int j = 0; j < N; ++j) {
    mu_sub(j) = mu(j * fraction);
  }
  return mu_sub;
}
  

} /* namespace FisherKPP */
