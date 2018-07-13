#ifndef __TPMESHBUILDER_H__
#define __TPMESHBUILDER_H__

#include <lf/base/static_vars.h>
#include "mesh_factory.h"

namespace lf::mesh::hybrid2d {

/**
 * @ brief Builder interface for creating structured meshes of rectangular
 * domains
 *
 * The design of the class compliees with the builder pattern.
 *
 * The geometry of the rectangular domain can be set by specifying the corners
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
      std::shared_ptr<mesh::MeshFactory> mesh_factory)
      : mesh_factory_(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimMesh() == 2,
        "TPTriagMeshBuilder can only construct meshes with DimMesh==2");
  }

  /**
   * @defgroup tpinit Initialization methods
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
  StructuredMeshBuilder &setNoXCells(size_type nxc) {
    no_of_x_cells_ = nxc;
    return *this;
  }
  StructuredMeshBuilder &setNoYCells(size_type nyc) {
    no_of_y_cells_ = nyc;
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
  /** mesh factory object that has to be supplied to the MeshBuilder */
  std::shared_ptr<mesh::MeshFactory> mesh_factory_;
  /** corners of rectangle defining the domain */
  Eigen::Vector2d bottom_left_corner_, top_right_corner_;
  /** Mesh resolution parameters */
  size_type no_of_x_cells_{0}, no_of_y_cells_{0};

};  // end class definition StructuredMeshBuilder

/**
 * @brief Implements a `MeshBuilder` that generates a triangular structured mesh
 *
 * Mesh generator for a triangular tensor product mesh covering a rectangle.
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
class TPTriagMeshBuilder : public StructuredMeshBuilder {
 public:
  /**
   * @brief Constructor: set factory object to be used by the builder
   *
   */
  explicit TPTriagMeshBuilder(std::shared_ptr<mesh::MeshFactory> mesh_factory)
      : StructuredMeshBuilder(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimWorld() == 2,
        "TPTriagMeshBuilder can only construct meshes with DimWorld==2");
  }

  /**
   * @brief actual construction of the mesh
   *
   */
  virtual std::shared_ptr<mesh::Mesh> Build() override;

 private:
  /**
   * @brief vertex index from grid position
   */
  size_type VertexIndex(size_type i, size_type j) const {
    return (i + j * (no_of_x_cells_ + 1));
  }

 public:
  /** Diagnostics control variable */
  static int output_ctrl_;
};  // end class definition TPTriagMeshBuilder

/**
 * @brief Implements a Builder for a tensor product grid (with rectangular
 * cells)
 *
 * The generated grid will be uniform, which means that all cells are congruent.
 * The geoemtry of the rectangular domain and the number of cells can be
 * specified by setting Builder state parameters.
 *
 */
class TPQuadMeshBuilder : public StructuredMeshBuilder {
 public:
  /**
   * @brief Constructor: set factory object to be used by the builder
   *
   */
  explicit TPQuadMeshBuilder(std::shared_ptr<mesh::MeshFactory> mesh_factory)
      : StructuredMeshBuilder(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimWorld() == 2,
        "TPQuadMeshBuilder can only construct meshes with DimWorld==2");
  }

  /**
   * @brief actual construction of the mesh
   *
   */
  virtual std::shared_ptr<mesh::Mesh> Build() override;

 private:
  /**
   * @brief vertex index from grid position
   */
  size_type VertexIndex(size_type i, size_type j) const {
    return (i + j * (no_of_x_cells_ + 1));
  }

 public:
  /** Diagnostics control variable */
  static int output_ctrl_;
};  // end class definition TPQuadMeshBuilder

}  // namespace lf::mesh::hybrid2d

#endif  //  __TPMESHBUILDER_H__
