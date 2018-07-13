/**
 * @file
 * @brief Declares the class Refinement Pattern
 * @author Raffael Casagrande
 * @date   2018-07-13 09:04:01
 * @copyright MIT License
 */

#ifndef __3178e8d1e7bf4366bcb00cdb4ebbf5fb
#define __3178e8d1e7bf4366bcb00cdb4ebbf5fb
#include "lf/base/base.h"

namespace lf::geometry {
/**
 * @brief Abstract interface class for topological local refinement
 *
 * This class defines a local topological refinement pattern by convex
 * polygon on an _integer lattice_ covering the reference element.
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

  RefinementPattern(lf::base::RefEl ref_el, lf::base::size_type lattice_const)
      : ref_el_(std::move(ref_el)), lattice_const_(lattice_const) {
    // Lattice constant N should be a multiple of six in order to be able to
    // define both the barycenter (lattice coordinates [N/3,N/3,N/3], and
    // the midpoints of edges. (lattice coordinates e.g. [N/2,N/2,0]).
    // Of course all lattice coordinates must be integers and must add up to N.
    LF_VERIFY_MSG(lattice_const % 6 == 0,
                  "Lattice constant should be multiple of 6");
  }

  lf::base::RefEl RefEl() const { return ref_el_; }
  lf::base::size_type LatticeConst() const { return lattice_const_; }
  /**
   * @brief provide number of child cells to be created by refinement
   */
  virtual lf::base::size_type noChildren() const = 0;
  /**
   * @brief provide lattice reference coordinates of vertices of child polygons
   *
   * the shaoe of the children of a cell is defined through a convex lattice
   * polygon in the reference element.
   *
   * @return vector of integer 2xn matrices containing the lattice coordinates
   * of the verticess of the child polygons. The length of this vector must
   * agree with the value returned by `noChildren()`. The integer entries of the
   * matrices must be non-negative and the column sums must be <= the lattice
   * constant.
   */
  virtual std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
  ChildPolygons() const = 0;

  /** Virtual Destructor */
  virtual ~RefinementPattern() = default;

 protected:
  lf::base::RefEl ref_el_;            /**< cell type */
  lf::base::size_type lattice_const_; /**< defines spacing of integer lattice */
};
}  // namespace lf::geometry

#endif  // __3178e8d1e7bf4366bcb00cdb4ebbf5fb
