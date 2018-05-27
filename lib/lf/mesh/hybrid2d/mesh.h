#ifndef __b43a044f002e459b8060206d5193d7ec
#define __b43a044f002e459b8060206d5193d7ec

#include <tuple>

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2d {

template <char CODIM>
class Entity;

class Mesh : public mesh::Mesh {
 public:
  /**
   * @brief Create a new instance of this mesh by directly specifying the
   *        nodes and elements. Note that it is preferred to use the
   *        MeshBuilder interface instead.
   * @param dim_world The dimension of the euclidean space in which the mesh
   *                  is embedded.
   * @param nodes     The coordinates of the nodes.
   * @param elements  A vector that describes the elements:
   *
   *
   * - `std::get<0>(elements[i])` contains the zero-based indices of the nodes
   *   that make up this element.
   * - `std::get<1>(elements[i])` contains the geometry mapping for that element
   *
   * @note This method will deduce the geometries of the segments and nodes
   *       from the element Geometry objects.
   */
  Mesh(
      char dim_world, std::vector<Eigen::VectorXd> nodes,
      std::vector<std::tuple<std::vector<size_type>, std::unique_ptr<geometry::Geometry>>>
          elements);

  char DimMesh() const override { return 2; }

  char DimWorld() const override { return dim_world_; }

  base::ForwardRange<const mesh::Entity> Entities(char codim) const override;

  size_type Size(char codim) const override {
    LF_ASSERT_MSG(codim >= 0, "codim negative.");
    LF_ASSERT_MSG(codim <= dim_world_, "codim > dimWorld.");

    switch (codim) {
      case 0:
        return entities0_.size();
      case 1:
        return entities1_.size();
      case 2:
        return entities2_.size();
      default:
        LF_VERIFY_MSG(false, "Something is horribyl wrong, codim = " +
                                 std::to_string(codim) + " is out of bounds.");
    }
  }

  size_type Index(const mesh::Entity& e) const override;

 private:
  char dim_world_;
  std::vector<Entity<0>> entities0_;
  std::vector<Entity<1>> entities1_;
  std::vector<Entity<2>> entities2_;

  template <char CODIM>
  friend class Entity;
};

}  // namespace lf::mesh::hybrid2d

#endif  // __b43a044f002e459b8060206d5193d7ec
