#ifndef THESIS_MESH_SHPERE_TRIAG_MESH_BUILDER_H
#define THESIS_MESH_SHPERE_TRIAG_MESH_BUILDER_H

/**
 * @file sphere_triag_mesh_builder.h
 * @brief Constructs a mesh on an annulu
 */

#define _USE_MATH_DEFINES

#include <lf/mesh/mesh.h>

namespace projects::hldo_sphere {

namespace mesh {

/**
 * @brief A mesh builder for regular 3-Sphere
 *
 * The diameter of the hole is controlled by the value passed to
 * setRadius().
 * The refinement level can be set with setRefinementLevel().
 * The refinement level indicates by how many times the edges of the
 * inital structure (octahedron, see https://de.wikipedia.org/wiki/Octahedron)
 * will be regualarly split. Note that this is only conceptual, the
 * implementation does not split egdes but directly build the desired sphere.
 * The radius defaults to 1.0 while refinement level defaults to 0.
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
   * The refinement level indicates by how many times the edges of the
   * inital structure (icosahedron, see
   * https://en.wikipedia.org/wiki/Icosahedron) will be split regualarly. Note
   * that this is only conceptual, the implementation does not split egdes but
   * directly build the desired sphere.
   */
  void setRefinementLevel(lf::base::size_type n) { refinement_level_ = n; }

  /**
   * @brief Build the mesh
   * @returns A shared pointer to a mesh in the form of a sphere. The radius is
   * set with `setRadius(r)` and the refinement level is set with
   * `setRefinementLevel(n)`.
   */
  std::shared_ptr<lf::mesh::Mesh> Build();

 private:
  double radius_;
  lf::base::size_type refinement_level_;
  std::unique_ptr<lf::mesh::MeshFactory> mesh_factory_;
};

}  // end namespace mesh

}  // namespace projects::hldo_sphere

#endif  // THESIS_MESH_SPHERE_TRIAG_MESH_BUILDER_H
