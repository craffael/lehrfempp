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

using size_type = lf::base::size_type;
using dim_t = lf::base::dim_t;
using sub_idx_t = lf::base::sub_idx_t;
using glb_idx_t = lf::base::glb_idx_t;
const unsigned int idx_nil = lf::base::kIdxNil;

class MeshFactory;

/** @brief Basis 2D mesh type compliant with abstract mesh interface
 *
 */
class Mesh : public mesh::Mesh {
 public:
  char DimMesh() const override { return 2; }
  char DimWorld() const override { return dim_world_; }

  base::ForwardRange<const mesh::Entity> Entities(char codim) const override;
  size_type Size(char codim) const override;
  size_type Index(const Entity& e) const override;
  const mesh::Entity* EntityByIndex(dim_t codim,
                                    glb_idx_t index) const override;
  bool Contains(const mesh::Entity& e) const override;

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
  /** @brief Auxliary array of cell (co-dim ==0 entities) pointers
   *
   * This array serves two purposes. It facilitates the construction
   * of a range covering all the cells. It is also required to retrieving
   * the entity of a specific codimension belonging to an index.
   */
  std::array<std::vector<const mesh::Entity*>, 3> entity_pointers_;

  /** @brief Data types for passing information about mesh intities */
  using GeometryPtr = std::unique_ptr<geometry::Geometry>;
  using NodeCoordList = std::vector<GeometryPtr>;
  using EdgeList =
      std::vector<std::pair<std::array<size_type, 2>, GeometryPtr>>;
  using CellList =
      std::vector<std::pair<std::array<size_type, 4>, GeometryPtr>>;

  /**
   * @brief Construction of mesh from information gathered in a MeshFactory
   * @param dim_world Dimension of the ambient space.
   * @param nodes sequential container of node coordinates
   * @param edges sequential container of pairs of
                  (i) vectors of indices of the nodes of an edge
                  (ii) pointers to the geometry object describing an edge
   * @param cells sequential container of pairs of
                  (i) vectors of indices of the nodes of a cell
                  (ii) pointers to the geometry object for the cell
   *
   * ### Shape guessing
   *
   * Shape information usually supplied by a unique pointer to a
   lf::geometry::Geometry
   * object attached to an entity may be missing for some of the entities. In
   this
   * case the constructor tries to reconstruct the shape from that of other
   entities.
   *
   * The main example is missing geometry information for nodes or edges. In
   this case
   * a call to the lf::geometry::SubGeometry() method of cell entities is used
   to
   * generate shape information for its lower-dimensional sub-entities. Note
   that
   * any adjacent cell can be used and no selection rule is given.
   *
   * Even cells without shape information may be passed. In this case the
   current implementation
   * builds a cell with straight edges of type lf::geometry::TriaO1 (in the case
   of a topological
   * triangle) or lf::geometry::QuadO1 (in the case of a quadrilateral) from the
   node positions.
   * The shape of edges is not taken into account.
   *
   * An extreme situation is marked by passing node positions as the only
   geometric information
   * together with topological node-cell incidence relationships. In this case
   the constructor
   * will build a mesh with straight edges throughout, type
   lf::geometry::SegmentO1.
   *
   * ### Missing entities
   *
   * @note not all edges of the mesh have to be passed in the `edges` argument;
   *       missing edges are inserted based on the information about the cells.
   *
   * @note the position of node information the `nodes` array and of cell
   *        information in the `cells` array, respectively,
   *        determines the interpretation of the index numbers,
   *        that is the n-th node in the container has index n-1.
   *
   */
  Mesh(dim_t dim_world, NodeCoordList nodes, EdgeList edges, CellList cells);

  friend class MeshFactory;

 public:
  /** @brief Diagnostics control variable */
  static int output_ctrl_;
};

}  // namespace lf::mesh::hybrid2dp

#endif  // __62731052ee4a4a2d9f256c2caac43835
