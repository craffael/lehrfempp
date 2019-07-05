#ifndef _LF_REFINEMENT_PAT_H_
#define _LF_REFINEMENT_PAT_H_

/**
 * @file refinement_pattern.h
 * @brief Definitions of data structures related to refinement of the mesh
 *
 */

#include <lf/mesh/mesh.h>

namespace lf::refinement {
// Convenient types inherited from base namespace
using lf::base::dim_t;
using lf::base::glb_idx_t;
using lf::base::size_type;
using lf::base::sub_idx_t;
const unsigned int idx_nil = lf::base::kIdxNil;

/**
 * @brief (possibly incomplete) list of refinement patterns
 *
 * _Non-symmetric_ refinement patterns have to be complemented by a subentity
 * index of an edge,  here called the **anchor** of the refinement.
 *
 * The only symmetric refinement patterns are regular refinement and barycentric
 * refinement.
 *
 * A refinement pattern has to be passed as first argument to the member
 * function ChildGeometry of a geometry object.
 *
 */
enum RefPat : int {
  rp_nil,          /**< no refinement */
  rp_copy,         /**< just copy the geometry of the current entity */
  rp_split,        /**< for edges: split into two halves */
  rp_bisect,       /**< single bisection, splitting refinement edge */
  rp_trisect,      /**< double bisection */
  rp_trisect_left, /**< double bisection, other triangle in second step */
  rp_quadsect,     /**< triple bisection, starting at edge 0 */
  rp_threeedge,    /**< split quad into one half and three triangles */
  rp_regular,      /**< standard "full" (regular) refinement */
  rp_barycentric   /**< barycentric refinement */
};

std::ostream& operator<<(std::ostream& o, const RefPat& refpat);

/**
 * @brief Class containing information about the refinement of a cell
 *
 * For explanations see xournal notes in Refinement.xoj
 */
class Hybrid2DRefinementPattern : public geometry::RefinementPattern {
  /** Constructors */
  /** @{ */
  Hybrid2DRefinementPattern(const Hybrid2DRefinementPattern&) = default;
  Hybrid2DRefinementPattern(Hybrid2DRefinementPattern&&) = default;
  Hybrid2DRefinementPattern& operator=(const Hybrid2DRefinementPattern&) =
      default;
  Hybrid2DRefinementPattern& operator=(Hybrid2DRefinementPattern&&) = default;
  /** @} */

  // TODO(ralfh):
  // We do not need the achor_set_ member, if we follow the convention that
  // lf::base::kIdxNil represents an invalid (= non-set) anchor index.

 public:
  /** @brief constructor */
  explicit Hybrid2DRefinementPattern(lf::base::RefEl ref_el)
      : geometry::RefinementPattern(ref_el),
        anchor_(idx_nil),
        ref_pat_(rp_nil),
        anchor_set_(false) {}

  /** @brief constructor */
  Hybrid2DRefinementPattern(lf::base::RefEl ref_el, RefPat ref_pat)
      : geometry::RefinementPattern(ref_el),
        anchor_(idx_nil),
        ref_pat_(ref_pat),
        anchor_set_(false) {}

  /** @brief constructor */
  Hybrid2DRefinementPattern(lf::base::RefEl ref_el, RefPat ref_pat,
                            lf::base::sub_idx_t anchor)
      : geometry::RefinementPattern(ref_el),
        anchor_(anchor),
        ref_pat_(ref_pat),
        anchor_set_(anchor != idx_nil) {
    if (anchor_set_) {
      LF_VERIFY_MSG(
          anchor < ref_el_.NumSubEntities(1),
          "Anchor " << anchor << " invalid for " << ref_el_.ToString());
    }
  }

  /**
   * @copydoc lf::geometry::RefinementPattern::noChildren()
   *
   * For a point: 0 (`rp_nil`), 1(`rp_copy`)
   */
  lf::base::size_type noChildren(lf::base::dim_t codim) const override;
  /**
   * @copydoc lf::geometry::RefinementPattern::ChildPolygons()
   *
   * ### Case of a point entity (dimension 0)
   *
   * The method always returns an single "matrix" of size 0xP, where
   * P=1 for `rp_copy`, and P=0 for `rp_nil`.
   *
   * ### Case of a segment entity (dimension 1)
   *
   * - codim = 0, request information about child segments.
   *              The method returns a Q-vector of 1xP integer matrices:
   *     + Q=0 for `rp_nil`
   *     + Q=1 for `rp_copy`, P=2, returns `[0 N]'
   *     + Q=2 for `rp_split`, P=2, returns '{[0 N/2],[N/2 N]}`
   *
   * - codim = 1: provide information about new _interior_ points created
   *              during refinement.
   * The method returns an empty vector unless `rp_split` is the refinement
   * pattern, In this case a 1x1 matrix `[N/2]` is returned.
   *
   * Here, `N` denotes the lattice constant, corresponding to the float
   * value 1.0 in reference coordinates
   *
   * ### Case of triangular cell entity (dimension 2)
   *
   * The method returns a Q-vector of 2xP integer matrices.
   *
   * - codim=0: request information about child cells
   *
   * Below the output is largely visualized by pictures. The big pink numbers
   * give the local index of the child cells, which is also the index in the
   * return vector. The small orange number indicate the local vertex indices in
   * the child cell.
   *
   *     + Q=0 for `rp_nil` (no refinement at all)
   *     + Q=1 for `rp_copy`, returns `{[0 N 0;0 0 N]}`
   *     + Q=2 for `rp_bisect`, anchor=0
   * @image html refinement_tria/rp_bisect_0.png width=400px
   *     + Q=2 for `rp_bisect`, anchor=1
   * @image html refinement_tria/rp_bisect_1.png width=400px
   *     + Q=2 for `rp_bisect`, anchor=2
   * @image html refinement_tria/rp_bisect_2.png width=400px
   *     + Q=3 for `rp_trisect`, anchor=0
   * @image html refinement_tria/rp_trisect_0_1.png width=400px
   *     + Q=3 for `rp_trisect`, anchor=1
   * @image html refinement_tria/rp_trisect_1_2.png width=400px
   *     + Q=3 for `rp_trisect`, anchor=2
   * @image html refinement_tria/rp_trisect_2_0.png width=400px
   *     + Q=3 for `rp_trisect_left`, anchor=0
   * @image html refinement_tria/rp_trisect_0_2.png width=400px
   *     + Q=3 for `rp_trisect_left`, anchor=1
   * @image html refinement_tria/rp_trisect_1_0.png width=400px
   *     + Q=3 for `rp_trisect_left`, anchor=2
   * @image html refinement_tria/rp_trisect_2_1.png width=400px
   *     + Q=4 for `rp_quadsect`, anchor=0
   * @image html refinement_tria/rp_quadsect_0.png width=400px
   *     + Q=4 for `rp_quadsect`, anchor=1
   * @image html refinement_tria/rp_quadsect_1.png width=400px
   *     + Q=4 for `rp_quadsect`, anchor=2
   * @image html refinement_tria/rp_quadsect_2.png width=400px
   *     + Q=4 for `rp_regular`:
   * @image html refinement_tria/rp_regular.png width=400px
   *     + Q=6 for `rp_barycentric`
   * @image html refinement_tria/rp_barycentric.png width=400px
   *
   * - codim=1: return information about (new) _interior_ child edges
   *
   * - codim=0: tell about new _interior_ child nodes
   *   Only for the refinement pattern `rp_barycentric` a new interior
   *   node is created (Q=1), Q=0 in all other cases.
   *
   * ### Case of quadrilateral cell (dimension 2)
   *
   * A Q-vector of 2xP integer matrices with lattice coordinates in their
   * columns is returned.
   *
   * - codim=0: tell about all child cells
   *
   * Below we visualize the output with large pink numbers indicating the
   * the index of a child cell in the return vector.
   *
   *    + Q=0 for `rp_nil` (no refinement)
   *    + Q=1 for `rp_copy` (quadrilateral is duplicated)
   *    + Q=2 for `rp_split'/`rp_bisect`, only visualized for anchor=0, the
   * other numberings arise from cyclic permutations.
   * @image html refinement_quad/rp_split_0.png width=500px
   *    + Q=3 for `rp_trisect`, image for anchor=0
   * @image html refinement_quad/rp_trisect_0.png width=500px
   *    + Q=4 for `rp_quadsect`, image for anchor=0
   * @image html refinement_quad/rp_quadsect_0.png width=500px
   *    + Q=4 for `rp_threeedge`
   *    + Q=4 for `rp_regular`/`rp_barycentric`
   *
   * - codim=1: obtain information about _interior_ child edges
   *
   * - codim=2: information about _interior_ child nodes
   *
   */
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> ChildPolygons(
      lf::base::dim_t codim) const override;

  /** @brief set local number of anchor edge */
  Hybrid2DRefinementPattern& setAnchor(lf::base::sub_idx_t anchor) {
    anchor_ = anchor;
    anchor_set_ = anchor != idx_nil;
    if (anchor_set_) {
      LF_VERIFY_MSG(
          anchor < ref_el_.NumSubEntities(1),
          "Anchor " << anchor << " invalid for " << ref_el_.ToString());
    }
    return (*this);
  }

  /** @brief set refinement pattern
   *
   * @note invalidates the local number of the anchor edge
   *
   * the calling sequence should be
   * `setRefPattern(..).setAnchor(..)`
   */
  Hybrid2DRefinementPattern& setRefPattern(RefPat ref_pat) {
    ref_pat_ = ref_pat;
    anchor_set_ = false;
    return *this;
  }

  /** @name Access methods
   * @{
   */
  lf::base::sub_idx_t anchor() const { return anchor_; }
  RefPat refpat() const {
    LF_VERIFY_MSG(
        !(((ref_pat_ == rp_bisect) || (ref_pat_ == rp_trisect) ||
           (ref_pat_ == rp_trisect_left) || (ref_pat_ == rp_quadsect) ||
           (ref_pat_ == rp_quadsect) || (ref_pat_ == rp_threeedge)) &&
          !anchor_set_),
        "ref pattern " << ref_pat_ << " needs anchor!");
    return ref_pat_;
  }
  /** @} */

  ~Hybrid2DRefinementPattern() override = default;

 private:
  lf::base::sub_idx_t anchor_; /**< local number of anchor edge */
  RefPat ref_pat_;             /**< refinement pattern */
  bool anchor_set_;            /**< flag indicating valid anchor */
};

}  // namespace lf::refinement

#endif
