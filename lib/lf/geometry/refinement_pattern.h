/**
 * @file
 * @brief Declares the class Refinement Pattern
 * @author Raffael Casagrande
 * @date   2018-07-13 09:04:01
 * @copyright MIT License
 */

#ifndef __3178e8d1e7bf4366bcb00cdb4ebbf5fb
#define __3178e8d1e7bf4366bcb00cdb4ebbf5fb
#include <utility>
#include "lf/base/base.h"

namespace lf::geometry {

/**
 * @brief Abstract interface class for encoding topological local refinement
 *
 * This class defines a local topological refinement pattern by convex
 * polygons on an _integer lattice_ covering the reference element.
 *
 * The rationale for using an integer lattice is the possibility of exact
 * arithmetic. It also emphasizes the hybrid character of refinement in
 * between topological and geometric operations.
 *
 */
class RefinementPattern {
 protected:
  RefinementPattern(const RefinementPattern&) = default;
  RefinementPattern(RefinementPattern&&) = default;
  RefinementPattern& operator=(const RefinementPattern&) = default;
  RefinementPattern& operator=(RefinementPattern&&) = default;

 public:
  explicit RefinementPattern(lf::base::RefEl ref_el)
      : ref_el_(std::move(ref_el)), lattice_const_(6) {}

  /** @brief Constructor fixing reference element and refinement resolution
   *
   * @param ref_el topological reference element
   */

  RefinementPattern(lf::base::RefEl ref_el, lf::base::size_type lattice_const)
      : ref_el_(std::move(ref_el)), lattice_const_(lattice_const) {
    // Lattice constant N should be a multiple of six in order to be able to
    // define both the barycenter (lattice coordinates [N/3,N/3,N/3], and
    // the midpoints of edges. (lattice coordinates e.g. [N/2,N/2,0]).
    // Of course all lattice coordinates must be integers and must add up to N.
    LF_VERIFY_MSG(lattice_const % 6 == 0,
                  "Lattice constant should be multiple of 6");
  }

  [[nodiscard]] lf::base::RefEl RefEl() const { return ref_el_; }
  [[nodiscard]] lf::base::size_type LatticeConst() const {
    return lattice_const_;
  }
  /**
   * @brief provide number of child entities of a given co-dimension
   * to be created by refinement
   */
  [[nodiscard]] virtual lf::base::size_type NumChildren(
      lf::base::dim_t codim) const = 0;
  /**
   * @brief provide lattice reference coordinates of vertices of child polygons
   *
   * @param codim _relative_ codimension of the children whose lattice polygons
   *        are requested.
   *
   * ### for a cell entity
   *
   * The shaoe of the children of relative co-dimension 0 of a cell is
   * defined through a convex lattice polygon in the reference element.
   * The children on co-dimension 1 are _interior_ edges. Their shape
   * is described by lattice segments. Children of relative co-dimension 2
   * are _interior_ points. Their position is given by a single lattice point.
   *
   * ### For a segment entity
   *
   * The shape of children with relative co-dimension 0 is given by
   * lattice intervals. Children with relative co-dimension 1 are
   * interior points and their location is given by single lattice points.
   *
   * @return vector of _integer matrices_ containing the **lattice coordinates**
   * of the verticess of the child polygons in their columns.
   * The size of the matrices is dxP, where d is the intrinsic dimension of the
   * entity, and P stands for the number of vertices of a particular child
   * entity.
   * The length of the returned vector must agree with the value returned
   * by `NumChildren()`
   * The integer entries of the matrices must be non-negative and the
   * column sums must be <= the lattice constant.
   */
  [[nodiscard]] virtual std::vector<
      Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
  ChildPolygons(lf::base::dim_t codim) const = 0;

  virtual ~RefinementPattern() = default;

 protected:
  // TODO(craffael): make this a pure virtual base class
  // NOLINTNEXTLINE
  lf::base::RefEl ref_el_; /**< cell type */
  // NOLINTNEXTLINE
  lf::base::size_type lattice_const_; /**< defines spacing of integer lattice */
};

/**
 * @brief Test whether a lattice polygon describes a logical parallelogram
 *
 * @param polygon an integer matrix whose column contain the lattice coordinates
 *        of the vertices of the polygon.
 *
 * A polygon passes the parallelogram test, if
 * -# it has four vertices
 * -# the vectors \f$x_1-x_0\f$ and \f$x_2-x_3\f$ agree.
 */
bool isParallelogram(
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>& polygon);

}  // namespace lf::geometry

#endif  // __3178e8d1e7bf4366bcb00cdb4ebbf5fb
