#ifndef __02a3dfa9ae3a4969b29d4c0ecfaa6ad9
#define __02a3dfa9ae3a4969b29d4c0ecfaa6ad9

/**
 * @file geometry.h
 * \brief Defines the Geometry interface and provides a number of
 *        classes that implement this interface + additional geometry related
 *        helper routines
 */
#include "geometry_interface.h"
#include "point.h"
#include "quad_o1.h"
#include "segment_o1.h"
#include "tria_o1.h"

namespace lf::geometry {

/**
 * @brief Abstract interface class for topological local refinement
 *
 * This class defines a local topological refinement pattern by convex
 * polygon on an _integer lattice_ covering the reference element.
 *
 */
class RefinementPattern {
 public:
  RefinementPattern(lf::base::RefEl ref_el)
      : ref_el_(ref_el), lattice_const_(6) {}

  RefinementPattern(lf::base::RefEl ref_el, lf::base::size_type lattice_const)
      : ref_el_(ref_el), lattice_const_(lattice_const) {
    // Lattice constant N should be a multiple of six in order to be able to
    // define both the barycenter (lattice coordinates [N/3,N/3,N/3], and
    // the midpoints of edges. (lattice coordinates e.g. [N/2,N/2,0]).
    // Of course all lattice coordinates must be integers and must add up to N.
    LF_VERIFY_MSG(lattice_const % 6 == 0,
                  "Lattice constant should be multiple of 6");
  }

  lf::base::RefEl RefEl(void) const { return ref_el_; }
  lf::base::size_type LatticeConst(void) const { return lattice_const_; }
  /**
   * @brief provide number of child cells to be created by refinement
   */
  virtual lf::base::size_type noChildren(void) const = 0;
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
  ChildPolygons(void) const = 0;

 protected:
  lf::base::RefEl ref_el_;            /**< cell type */
  lf::base::size_type lattice_const_; /**< defines spacing of integer lattice */
};

}  // namespace lf::geometry

#endif  // __02a3dfa9ae3a4969b29d4c0ecfaa6ad9
