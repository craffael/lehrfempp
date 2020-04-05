#ifndef THESIS_MESH_ANNULUS_TRIAG_MESH_BUILDER_H
#define THESIS_MESH_ANNULUS_TRIAG_MESH_BUILDER_H

/**
 * @file annulus_triag_mesh_builder.h
 * @brief Constructs a mesh on an annulus
 */

#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/structured_mesh_builder.h>

namespace projects::ipdg_stokes {

namespace mesh {

/**
 * @brief A mesh builder for disks with a hole in the middle
 *
 * The diameter of the hole is controlled by the value passed to
 * `setInnerRadius(r)` while the diameter of the disk is controlled by
 * the value passed to `setOutrRadius(r)`. The position of the center
 * can be controlled by setting `setCenterPoint(x, y)`.
 * The inner radius defaults to 0.25 while the outer radius defaults
 * to 1. The center is at the origin by default.
 */
class AnnulusTriagMeshBuilder {
 public:
  /**
   * @brief Constructor
   * @param mesh_factory A unique pointer to a mesh factory object used for the
   * assembly of the mesh
   */
  explicit AnnulusTriagMeshBuilder(
      std::unique_ptr<lf::mesh::MeshFactory> mesh_factory)
      : inner_radius_(0.25),
        outer_radius_(1),
        num_radial_cells_(1),
        num_angular_cells_(4),
        center_point_(Eigen::Vector2d::Zero()),
        mesh_factory_(std::move(mesh_factory)) {
    // Nothing to do here
  }

  /**
   * @breif Sets the radius of the hole
   * @param r The radius of the hole
   */
  void setInnerRadius(double r) { inner_radius_ = r; }

  /**
   * @brief Set the radius of the disk
   * @param r The radius of the disk
   */
  void setOuterRadius(double r) { outer_radius_ = r; }

  /**
   * @brief Set the location of the center
   * @param x The x coordinate of the new center of the annulus
   * @param y The y coordinate of the new center of the annulus
   */
  void setCenterPoint(double x, double y) { center_point_ << x, y; }

  /**
   * @brief Set the number of cells in the radial direction
   * @param n The number of cells in the radial direction
   */
  void setNumRadialCells(lf::base::size_type n) { num_radial_cells_ = n; }

  /**
   * @brief Set the number of cells in the angular direction
   * @param n The nummber of cells in the angular direction
   */
  void setNumAngularCells(lf::base::size_type n) { num_angular_cells_ = n; }

  /**
   * @brief Build the mesh
   * @returns A shared pointer to a mesh in the form of a disk with a hole in
   * the middle
   */
  std::shared_ptr<lf::mesh::Mesh> Build();

 private:
  double inner_radius_;
  double outer_radius_;
  lf::base::size_type num_radial_cells_;
  lf::base::size_type num_angular_cells_;
  Eigen::Vector2d center_point_;
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_;
};

}  // end namespace mesh

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_MESH_ANNULUS_TRIAG_MESH_BUILDER_H
