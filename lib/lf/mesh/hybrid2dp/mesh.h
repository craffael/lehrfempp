/**
 * @file
 * @brief Declares the Hybrid2dP Mesh class
 * @author Raffael Casagrande
 * @date   2018-06-22 03:47:33
 * @copyright MIT License
 */

#ifndef __62731052ee4a4a2d9f256c2caac43835
#define __62731052ee4a4a2d9f256c2caac43835

#include <lf/base/static_vars.h>
#include <lf/mesh/mesh.h>
#include "lf/mesh/utils/print_info.h"
#include "point.h"
#include "quad.h"
#include "segment.h"
#include "triangle.h"

namespace lf::mesh::hybrid2dp {

class MeshFactory;

/** @brief Basis 2D mesh type compliant with abstract mesh interface
 *
 */
class Mesh : public mesh::Mesh {
 public:
  using dim_t = lf::base::dim_t;
  using sub_idx_t = lf::base::sub_idx_t;
  using glb_idx_t = lf::base::glb_idx_t;
  char DimMesh() const override { return 2; }
  char DimWorld() const override { return dim_world_; }

  base::ForwardRange<const Entity> Entities(char codim) const override;
  size_type Size(char codim) const override;
  size_type Index(const Entity& e) const override;
  bool Contains(const Entity& e) const override;

 private:
  dim_t dim_world_{};
  /** @brief array of 0-dimensional entity object of co-dimension 2 */
  std::vector<hybrid2dp::Point> points_;
  /** @brief array of 1-dimensional entity object of co-dimension 1 */
  std::vector<hybrid2dp::Segment> segments_;
  /** @brief array of triangular cell objects, oo-dimension 0 */
  std::vector<hybrid2dp::Triangle> trias_;
  /** @brief array of quadrilateral cell objects, oo-dimension 0 */
  std::vector<hybrid2dp::Quadrilateral> quads_;
  /** @brief Auxliary array of cell pointers */
  std::vector<const mesh::Entity*> cell_pointers_;

  /** @brief Data types for passing information about mesh intities */
  using NodeCoordList = std::vector<Eigen::VectorXd>;
  using GeometryPtr = std::unique_ptr<geometry::Geometry>;
  using EdgeList =
      std::vector<std::pair<std::array<Mesh::size_type, 2>, GeometryPtr>>;
  using CellList =
      std::vector<std::pair<std::array<Mesh::size_type, 4>, GeometryPtr>>;

  /**
   * @brief Construction of mesh from information gathered in a MeshFactory
   * @param nodes sequential container of node coordinates
   * @param edges sequential container of pairs of
                  (i) vectors of indices of the nodes of an edge
                  (ii) pointers to the geometry object describing an edge
   * @param cells sequential container of pairs of
                  (i) vectors of indices of the nodes of a cell
                  (ii) pointers to the geometry object for the cell
   * @return TBD
   * @note the position of node information the `nodes` array and of cell
   *        information in the `cells` array, respectively,
   *        determines the interpretation of the index numbers,
   *        that is the n-th node in the container has index n-1.
   *
   */
  Mesh(dim_t dim_world, const NodeCoordList& nodes, EdgeList edges,
       CellList cells);

  friend class MeshFactory;

 public:
  /** @brief diagnostics control variable */
  static int output_ctrl_;
};

inline std::ostream& operator<<(std::ostream& stream, const Mesh& mesh){
    //stream << "mesh object";
    //utils::PrintInfo(mesh, stream);
}

}  // namespace lf::mesh::hybrid2dp

#endif  // __62731052ee4a4a2d9f256c2caac43835
