#ifndef __b43a044f002e459b8060206d5193d7ec
#define __b43a044f002e459b8060206d5193d7ec

#include <tuple>

#include <lf/mesh/mesh.h>
//#include "lf/mesh/utils/utils.h"

namespace lf::mesh::hybrid2d {

template <char CODIM>
class Entity;

class Mesh : public mesh::Mesh {
 public:
  char DimMesh() const override { return 2; }

  char DimWorld() const override { return dim_world_; }

  base::ForwardRange<const mesh::Entity> Entities(char codim) const override;

  size_type Size(char codim) const override;

  size_type Index(const mesh::Entity& e) const override;

  const mesh::Entity* EntityByIndex(dim_t codim,
                                    base::glb_idx_t index) const override;

  bool Contains(const mesh::Entity& e) const override;

 private:
  char dim_world_;
  std::vector<Entity<0>> entities0_;
  std::vector<Entity<1>> entities1_;
  std::vector<Entity<2>> entities2_;

  /**
   * @brief Create a new instance of this mesh by directly specifying the
   *        nodes and elements. Do not call this method directly, instead
   *        use a mesh::MeshBuilder to do this.
   * @param dim_world The dimension of the euclidean space in which the mesh
   *                  is embedded.
   * @param nodes     The coordinates of the nodes.
   * @param edges     Edges that have been added specifically by the user.
   * @param elements  A vector that describes the elements:
   *
   *
   * - `std::get<0>(elements[i])` contains the zero-based indices of the nodes
   *   that make up this element. If `std::get<0>(elements[i])[3]==-1`, the
   *   element is a tria.
   * - `std::get<1>(elements[i])` contains the geometry mapping for that element
   * - `std::get<0>(edges[i])` contains the two nodes that delimit the edge.
   * - `std::get<1>(edges[i])` contains the geometry mapping for the edge.
   *
   * @note This method will deduce the geometries of the segments and nodes
   *       from the element Geometry objects.
   */

  Mesh(char dim_world, std::vector<Eigen::VectorXd> nodes,
       std::vector<std::tuple<std::array<size_type, 2>,
                              std::unique_ptr<geometry::Geometry>>>
           edges,
       std::vector<std::tuple<std::array<size_type, 4>,
                              std::unique_ptr<geometry::Geometry>>>
           elements);

  template <char CODIM>
  friend class Entity;
  friend class MeshFactory;
};

}  // namespace lf::mesh::hybrid2d

#endif  // __b43a044f002e459b8060206d5193d7ec
