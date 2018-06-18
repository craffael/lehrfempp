#ifndef __MESHBUILDER_H__
#define __MESHBUILDER_H__

#include "mesh_factory.h"

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
class TPTriagMeshBuilder: public mesh::MeshBuilder {
public:
  /**
   * @brief Constructor: does nothing
   *
   */
  TPTriagMeshBuilder():no_of_x_cells_(0),no_of_y_cells_(0) {} 
   
  /** @copydoc Mesh::DimWorld */
  dim_t DimWorld() const override {return 2; }

  /** 
   * @brief 2D hybrid meshes are meant to model 2D manifolds
   *
   */
  dim_t DimMesh() const override {return 2; }

  /**
   * @brief Initialization methods
   * @{
   */
  template <typename VECTOR>
  TPTriagMeshBuilder &setBottomLeftCorner(const VECTOR &&blc) {
    bottom_left_corner_ << blc[0],blc[1]; return *this;
  }
  template <typename VECTOR>
  TPTriagMeshBuilder &setTopRightCorner(const VECTOR &&trc) {
    top_right_corner_ << trc[0],trc[1]; return *this;
  }
  TPTriagMeshBuilder &setNoXCells(size_type nxc) { 
    no_of_x_cells_ = nxc; return *this;
  }
  TPTriagMeshBuilder &setNoYCells(size_type nyc) { 
    no_of_y_cells_ = nyc; return *this;
  }
  /** @} */
  
  /** 
   * @brief actual construction of the mesh 
   *
   */
  std::unique_ptr<mesh::Mesh> Build() override;
private:
  /**
   * @brief vertex index from grid position
   */
  inline size_type VertexIndex(size_type i,size_type j) const {
    return (i+j*(no_of_x_cells_+1));
  }
  
private:
  Eigen::Vector2d bottom_left_corner_,top_right_corner_;
  size_type no_of_x_cells_,no_of_y_cells_;
};

#endif //  __MESHBUILDER_H__
