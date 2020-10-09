#ifndef STRUCTURED_MESH_BUILDER_H
#define STRUCTURED_MESH_BUILDER_H

#include "lf/mesh/mesh_factory.h"

namespace lf::mesh::utils {

/**
 * @ brief Builder interface for creating structured meshes on rectangular
 * domains
 *
 * The design of the class complies with the builder pattern.
 *
 * The geometry of the rectangular domain can be set by specifying the corners.
 */
class StructuredMeshBuilder {
 public:
  using size_type = mesh::Mesh::size_type;
  using dim_t = base::RefEl::dim_t;
  /*
   * @brief Constructor: set factory object to be used by the builder
   *
   */
  explicit StructuredMeshBuilder(
      std::unique_ptr<mesh::MeshFactory> mesh_factory)
      : mesh_factory_(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimMesh() == 2,
        "TPTriagMeshBuilder can only construct meshes with DimMesh==2");
  }

  /**
   * @name Initialization methods
   * @{
   */
  template <typename VECTOR>
  StructuredMeshBuilder &setBottomLeftCorner(const VECTOR &blc) {
    bottom_left_corner_ << blc[0], blc[1];
    return *this;
  }
  StructuredMeshBuilder &setBottomLeftCorner(double x0, double x1) {
    bottom_left_corner_ << x0, x1;
    return *this;
  }
  template <typename VECTOR>
  StructuredMeshBuilder &setTopRightCorner(const VECTOR &&trc) {
    top_right_corner_ << trc[0], trc[1];
    return *this;
  }
  StructuredMeshBuilder &setTopRightCorner(double x0, double x1) {
    top_right_corner_ << x0, x1;
    return *this;
  }
  StructuredMeshBuilder &setNumXCells(size_type nxc) {
    num_of_x_cells_ = nxc;
    return *this;
  }
  StructuredMeshBuilder &setNumYCells(size_type nyc) {
    num_of_y_cells_ = nyc;
    return *this;
  }
  /** @} */

  /**
   * @brief Interface for the actual construction of the mesh
   *
   * This method has to be implemented by derived classes
   *
   * @return a _shared_ pointer to the newly created mesh. This pointer can
   *         be copied freely and the mesh will exists until the last pointer to
   *         ceases to exist.
   */
  virtual std::shared_ptr<mesh::Mesh> Build() = 0;

 protected:
  // TODO(craffael): find a better way than implementation inheritance.
  /** mesh factory object that has to be supplied to the MeshBuilder */
  std::unique_ptr<mesh::MeshFactory> mesh_factory_;  // NOLINT
  /** corners of rectangle defining the domain */
  Eigen::Vector2d bottom_left_corner_, top_right_corner_;  // NOLINT
  /** Mesh resolution parameters */
  size_type num_of_x_cells_{0}, num_of_y_cells_{0};  // NOLINT
};  // end class definition StructuredMeshBuilder

}  // namespace lf::mesh::utils

#endif /* STRUCTURED_MESH_BUILDER_H */
