#ifndef __e98a803fac5b430a8ff634ceb2f809a1
#define __e98a803fac5b430a8ff634ceb2f809a1

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements mesh::MeshFactory interface and can be used to construct
 *        a hybrid mesh with `dimMesh=2`.
 *
 * A planar triangular mesh with affine triangles as cells can be completely
 * specified by giving the list of vertex coordinates and a list of triangles in
 * the form of a 3-tuple of vertex indices.
 */
class MeshFactory : public mesh::MeshFactory {
 public:
  /**
   * @brief Construct a new builder that can be used to construct a new hybrid2d
   *        mesh.
   * @param dim_world The dimension of the euclidean space in which the
   *                  mesh is embedded.
   */
  explicit MeshFactory(dim_t dim_world)
      : dim_world_(dim_world), built_(false) {}

  /** @copydoc Mesh::DimWorld() */
  dim_t DimWorld() const override { return dim_world_; }

  /** @copydoc Mesh::DimMesh() */
  dim_t DimMesh() const override { return 2; }

  /** @copydoc Mesh::AddPoint() */
  size_type AddPoint(coord_t coord) override;

  /** @copydoc Mesh::AddPoint() */
  size_type AddPoint(std::unique_ptr<geometry::Geometry>&& geometry) override;

  /** @copydoc Mesh::AddEntity() */
  size_type AddEntity(base::RefEl ref_el,
                      const base::ForwardRange<const size_type>& nodes,
                      std::unique_ptr<geometry::Geometry>&& geometry) override;

  /** @copydoc Mesh::Build() */
  std::shared_ptr<mesh::Mesh> Build() override;

 private:
  dim_t dim_world_;  // dimension of ambient space
  bool built_;
  std::vector<Eigen::VectorXd> nodes_;
  std::vector<
      std::tuple<std::array<size_type, 2>, std::unique_ptr<geometry::Geometry>>>
      edges_;
  std::vector<
      std::tuple<std::array<size_type, 4>, std::unique_ptr<geometry::Geometry>>>
      elements_;
};

}  // namespace lf::mesh::hybrid2d

#endif  // __e98a803fac5b430a8ff634ceb2f809a1
