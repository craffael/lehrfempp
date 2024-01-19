#ifndef HLDO_SPHERE_MESH_BUILDER_H
#define HLDO_SPHERE_MESH_BUILDER_H

/**
 * @file sphere_triag_mesh_builder.h
 * @brief Constructs a mesh on an annulu
 */

#define _USE_MATH_DEFINES

#include <lf/mesh/mesh.h>

namespace projects::hldo_sphere {

/**
 * @brief Contains a generator for triangluations of the surface of the 3-sphere
 */
namespace mesh {

/**
 * @brief A mesh builder for regular 3-Sphere
 *
 * The radius of the Sphere is controlled by setRadius().
 * The refinement level can be set with setRefinementLevel().
 * The refinement level indicates how fine the mesh should be
 * inital structure (octahedron, see https://de.wikipedia.org/wiki/Octahedron)
 * The radius defaults to 1.0 while refinement level defaults to 0.
 *
 * Details regarding the implementation and its concepts can be found in the
 * thesis `Hodge Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * section 4.1.
 */
class SphereTriagMeshBuilder {
 public:
  /**
   * @brief Constructor
   * @param mesh_factory A unique pointer to a mesh factory object used for the
   * assembly of the mesh
   */
  explicit SphereTriagMeshBuilder(
      std::unique_ptr<lf::mesh::MeshFactory> mesh_factory)
      : radius_(1.0),
        refinement_level_(0),
        mesh_factory_(std::move(mesh_factory)) {
    // Nothing to do here
  }

  /**
   * @brief Sets the radius
   * @param r The radius of the sphere
   *
   * requries r > 0
   */
  void setRadius(double r) {
    LF_ASSERT_MSG(0 < r, "radius needs to be positive");
    radius_ = r;
  }

  /**
   * @brief Set the refinement level
   * @param n The nummber of refinements
   *
   * The refinement level indicates by how fine the mesh should become.
   * Inital structure (octahedron, see
   * https://en.wikipedia.org/wiki/Octahedron) on which new vertices are
   * introduced for every refinement level.
   *
   * For refinement level l the mesh will have
   * @f$
   *   2^{2l + 3}
   * @f$
   * cells,
   * @f$
   *    3 \cdot 2^{2l + 2}
   * @f$
   * edges and
   * @f$
   *    2^{2l + 2} + 2
   * @f$
   * vertices.
   *
   */
  void setRefinementLevel(lf::base::size_type n) { refinement_level_ = n; }

  /**
   * @brief Build the mesh
   * @returns A shared pointer to a mesh in the form of a sphere. The radius is
   * set with `setRadius(r)` and the refinement level is set with
   * `setRefinementLevel(n)`.
   */
  std::shared_ptr<lf::mesh::Mesh> Build() const;

 private:
  double radius_;
  lf::base::size_type refinement_level_;
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_;
};

}  // end namespace mesh

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_MESH_BUILDER_H
