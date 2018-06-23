#ifndef __MESHBUILDER_H__
#define __MESHBUILDER_H__

#include "mesh_factory.h"

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements a `MeshBuilder` that generates a triangular structured mesh
 *
 * Mesh generator for a triangular tensor product mesh covering a rectangle.
 *
 * The design of the class compliees with the builder pattern.
 *
 * A triangular tensor product mesh of a rectangular domain is built by
 * subdividing the domain into equal squares and splitting each of them
 * along the diagonal.
 *
 * The nodes and cells of the resulting mesh are numbered lexikographically
 * from bottom left to  top right. The super-diagonal cells have even indices
 * the sub-diagonal cells odd indices.
 *
 * The horizontal edges are numbered first, the vertical edges next; both groups
 * lexikographically.
 */
class TPTriagMeshBuilder {
 public:
  using size_type = mesh::Mesh::size_type;
  using dim_t = base::RefEl::dim_t;
  /**
   * @brief Constructor: does nothing
   *
   */
  explicit TPTriagMeshBuilder(std::shared_ptr<MeshFactory> mesh_factory)
      : mesh_factory_(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimMesh() == 2,
        "TPTriagMeshBuilder can only construct meshes with DimMesh==2");
    LF_ASSERT_MSG(
        mesh_factory_->DimWorld() == 2,
        "TPTriagMeshBuilder can only construct meshes with DimWorld==2");
  }

  /**
   * @brief Initialization methods
   * @{
   */
  template <typename VECTOR>
  TPTriagMeshBuilder &setBottomLeftCorner(const VECTOR &blc) {
    bottom_left_corner_ << blc[0], blc[1];
    return *this;
  }
  TPTriagMeshBuilder &setBottomLeftCorner(double x0, double x1) {
    bottom_left_corner_ << x0, x1;
    return *this;
  }
  template <typename VECTOR>
  TPTriagMeshBuilder &setTopRightCorner(const VECTOR &&trc) {
    top_right_corner_ << trc[0], trc[1];
    return *this;
  }
  TPTriagMeshBuilder &setTopRightCorner(double x0, double x1) {
    top_right_corner_ << x0, x1;
    return *this;
  }
  TPTriagMeshBuilder &setNoXCells(size_type nxc) {
    no_of_x_cells_ = nxc;
    return *this;
  }
  TPTriagMeshBuilder &setNoYCells(size_type nyc) {
    no_of_y_cells_ = nyc;
    return *this;
  }
  /** @} */

  /**
   * @brief actual construction of the mesh
   *
   */
  std::shared_ptr<mesh::Mesh> Build();

 private:
  /**
   * @brief vertex index from grid position
   */
  size_type VertexIndex(size_type i, size_type j) const {
    return (i + j * (no_of_x_cells_ + 1));
  }

 private:
  std::shared_ptr<mesh::MeshFactory> mesh_factory_;
  Eigen::Vector2d bottom_left_corner_, top_right_corner_;
  size_type no_of_x_cells_{0}, no_of_y_cells_{0};
};

}  // namespace lf::mesh::hybrid2d

#endif  //  __MESHBUILDER_H__
