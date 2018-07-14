#ifndef _LF_REFINEMENT_H_
#define _LF_REFINEMENT_H_

/**
 * @file refinement.h
 * @brief Definitions of data structures related to refinement of the mesh
 *
 */

#include <lf/mesh/mesh.h>

namespace lf::refinement {
/**
 * @brief (possibly incomplete) list of refinement patterns
 *
 * _Non-symmetric_ refinement patterns have to be complemented by a subentity
 * index of an edge,  here called the anchor of the refinement.
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

 public:
  /** @brief constructor */
  explicit Hybrid2DRefinementPattern(lf::base::RefEl ref_el)
      : geometry::RefinementPattern(ref_el),
        anchor_(-1),
        ref_pat_(rp_nil),
        anchor_set_(false) {}

  /** @brief constructor */
  Hybrid2DRefinementPattern(lf::base::RefEl ref_el, RefPat ref_pat)
      : geometry::RefinementPattern(ref_el),
        anchor_(-1),
        ref_pat_(ref_pat),
        anchor_set_(false) {}

  /** @brief constructor */
  Hybrid2DRefinementPattern(lf::base::RefEl ref_el, RefPat ref_pat,
                            lf::base::sub_idx_t anchor)
      : geometry::RefinementPattern(ref_el),
        anchor_(anchor),
        ref_pat_(ref_pat),
        anchor_set_(true) {
    LF_VERIFY_MSG(anchor < ref_el_.NumSubEntities(1),
                  "Anchor " << anchor << " invalid for " << ref_el_.ToString());
  }

  /**
   * @copydoc RefinementPattern::noChildren
   */
  virtual lf::base::size_type noChildren(lf::base::dim_t codim) const;
  /**
   * @copydoc RefinementPattern::ChildPolygons
   */
  virtual std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
  ChildPolygons(lf::base::dim_t codim) const;

  /** @brief set local number of anchor edge */
  Hybrid2DRefinementPattern& setAnchor(lf::base::sub_idx_t anchor) {
    LF_VERIFY_MSG(anchor < ref_el_.NumSubEntities(1),
                  "Anchor " << anchor << " invalid for " << ref_el_.ToString());
    anchor_ = anchor;
    anchor_set_ = true;
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

  /** @defgroup getRefPat
   * @brief Access methods
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

  virtual ~Hybrid2DRefinementPattern() override = default;

 private:
  lf::base::sub_idx_t anchor_; /**< local number of anchor edge */
  RefPat ref_pat_;             /**< refinement pattern */
  bool anchor_set_;            /**< flag indicating valid anchor */
};

}  // namespace lf::refinement

#endif
