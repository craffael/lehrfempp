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

  /** @copydoc Mesh::DimWorld */
  dim_t DimWorld() const override { return dim_world_; }

  /**
   * @brief 2D hybrid meshes are meant to model 2D manifolds
   *
   */
  dim_t DimMesh() const override { return 2; }

  /**
   * @brief register the coordinates of another point
   * @param coord a dynamic Eigen vector of size DimWorld containing node
   * coordinates
   *
   */
  size_type AddPoint(coord_t coord) override;

  /*
   * @brief register another triangle
   * @param nodes a 3-tuple of node *indices*
   * @param geometry a description of the shape of a cell
   *
   */
  size_type AddElement(const base::ForwardRange<const size_type>& nodes,
                       std::unique_ptr<geometry::Geometry>&& geometry) override;

  /**
   * @brief actual construction of the mesh
   *
   */
  std::unique_ptr<mesh::Mesh> Build() override;

 private:
  dim_t dim_world_;  // dimension of ambient space
  bool built_;
  std::vector<Eigen::VectorXd> nodes_;
  std::vector<
      std::tuple<std::vector<size_type>, std::unique_ptr<geometry::Geometry>>>
      elements_;
};

}  // namespace lf::mesh::hybrid2d

#endif  // __e98a803fac5b430a8ff634ceb2f809a1
