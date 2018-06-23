/**
 * @file
 * @brief Declares the Hybrid2dP Mesh class
 * @author Raffael Casagrande
 * @date   2018-06-22 03:47:33
 * @copyright MIT License
 */

#ifndef __62731052ee4a4a2d9f256c2caac43835
#define __62731052ee4a4a2d9f256c2caac43835

#include <lf/mesh/mesh.h>
#include "point.h"
#include "quad.h"
#include "segment.h"
#include "triangle.h"

namespace lf::mesh::hybrid2dp {

class MeshFactory;

class Mesh : public mesh::Mesh {
 public:
  using dim_t = base::RefEl::dim_t;
  char DimMesh() const override { return 2; }
  char DimWorld() const override { return dim_world_; }

  base::ForwardRange<const Entity> Entities(char codim) const override;
  size_type Size(char codim) const override;
  size_type Index(const Entity& e) const override;

 private:
  dim_t dim_world_{};
  std::vector<Point> points_;
  std::vector<Segment> segments_;
  std::vector<Triangle> trias_;
  std::vector<Quadrilateral> quads_;

  // Data types for passing information about mesh intities
  using GeometryPtr = std::unique_ptr<geometry::Geometry>;
  using EdgeList =
      std::vector<std::pair<std::array<Mesh::size_type, 2>, GeometryPtr>>;
  using CellList =
      std::vector<std::pair<std::vector<Mesh::size_type>, GeometryPtr>>;

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
  Mesh(dim_t dim_world,
       std::vector<Eigen::VectorXd> nodes,
       EdgeList edges,
       CellList cells);
  
  friend class MeshFactory;
};

}  // namespace lf::mesh::hybrid2dp

#endif  // __62731052ee4a4a2d9f256c2caac43835
