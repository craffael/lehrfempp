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
  class LocalRefinePattern {
  public:
    LocalRefinePattern(lf::base::RefEl ref_el):
      ref_el_(ref_el),lattice_const_(6) {}
    
    LocalRefinePattern(lf::base::RefEl ref_el,lf::base::size_type lattice_const):
      ref_el_(ref_el),lattice_const_(lattice_const) {
      LF_VERIFY_MSG(lattice_const %6 == 0,"Lattice constant should be multiple of 6");
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
     * the shaoe of the children of a cell is defined through a convex lattice polygon
     * in the reference element. 
     *
     * @return vector of integer 2xn matrices containing the lattice coordinates
     * of the verticess of the child polygons. The length of this vector must agree
     * with the value returned by `noChildren()`. The integer entries of the matrices
     * must be non-negative and the column sums must be <= the lattice constant.
     */
    virtual std::vector<Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>>
    ChildPolygons(void) const = 0;
  protected:
    lf::base::RefEl ref_el_;     /**< cell type */
    lf::base::size_type lattice_const_; /**< defines spacing of integer lattice */
  };
  
 /**
 * @brief (possibly incomplete) list of refinement patterns
 *
 * _Non-symmetric_ refinement patterns have to be complemented by a subentity
 * index of an edge,  here called the anchor of the refinement.
 *
 * The only symmetric refinement patterns are regular refinement and barycentric 
 * refinement. 
 *
 * A refinement pattern has to be passed as first argument to the member function
 * ChildGeometry of a geometry object.
 *
 */
enum RefPat: int {
  rp_nil,      /**< no refinement */
  rp_copy,     /**< just copy the geometry of the current entity */
  rp_split,    /**< for edges: split into two halves */
  rp_bisect,    /**< single bisection, splitting refinement edge */
  rp_trisect,  /**< double bisection */
  rp_trisect_left, /**< double bisection, other triangle in second step */
  rp_quadsect,   /**< triple bisection, starting at edge 0 */
  rp_threeedge,  /**< split quad into one half and three triangles */
  rp_regular,    /**< standard "full" (regular) refinement */
  rp_barycentric /**< barycentric refinement */  
};

/** 
 * @brief Class containing information about the refinement of a cell
 * 
 * For explanations see xournal notes in Refinement.xoj
 */ 
  class RefinementPattern: public LocalRefinePattern {
public:
  /** @brief constructor */
  explicit  RefinementPattern(lf::base::RefEl ref_el):
    LocalRefinePattern(ref_el),anchor_(-1),ref_pat_(rp_nil),anchor_set_(false) {}
  
  /** @brief constructor */
  RefinementPattern(lf::base::RefEl ref_el,RefPat ref_pat):
    LocalRefinePattern(ref_el),anchor_(-1),ref_pat_(ref_pat),anchor_set_(false) {}
  
  /** @brief constructor */
  RefinementPattern(lf::base::RefEl ref_el,RefPat ref_pat,lf::base::sub_idx_t anchor):
    LocalRefinePattern(ref_el),anchor_(anchor),ref_pat_(ref_pat),anchor_set_(true)
  {
    LF_VERIFY_MSG(anchor < ref_el_.NumSubEntities(1),
		  "Anchor " << anchor << " invalid for " << ref_el_.ToString());
  }

    /**
     * @copydoc LocalRefinePattern::noChildren 
     */
    virtual lf::base::size_type noChildren(void) const;
    /**
     * @copydoc LocalRefinePattern::ChildPolygons
     */
    virtual std::vector<Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic>>
    ChildPolygons(void) const;

		       
  /** @brief set local number of anchor edge */
  RefinementPattern &setAnchor(lf::base::sub_idx_t anchor) {
    LF_VERIFY_MSG(anchor < ref_el_.NumSubEntities(1),
		  "Anchor " << anchor << " invalid for " << ref_el_.ToString());
    anchor_ = anchor; anchor_set_ = true;
    return (*this);
  }

  /** @brief set refinement pattern 
   *
   * @note invalidates the local number of the anchor edge
   *
   * the calling sequence should be
   * `setRefPattern(..).setAnchor(..)` 
   */
  RefinementPattern &setRefPattern(RefPat ref_pat) {
    ref_pat_ = ref_pat; anchor_set_ = false; }

  /** @defgroup getRefPat
   * @brief Access methods
   * @{
   */
  lf::base::sub_idx_t anchor(void) const { return anchor_; }
  RefPat refpat(void) const {
    LF_VERIFY_MSG(!(((ref_pat_ == rp_bisect) ||
		     (ref_pat_ == rp_trisect) ||
		     (ref_pat_ == rp_trisect_left) ||
		     (ref_pat_ == rp_quadsect) ||
		     (ref_pat_ == rp_quadsect) ||
		     (ref_pat_ == rp_threeedge)) && !anchor_set_),
		  "ref pattern " << (int)ref_pat_ << " needs anchor!");
    return ref_pat_; }
  /** @} */

  virtual ~RefinementPattern(void) = default;
private:
  lf::base::sub_idx_t anchor_; /**< local number of anchor edge */
  RefPat ref_pat_;             /**< refinement pattern */
  bool anchor_set_; /**< flag indicating valid anchor */
};

} // end namespace

#endif  // __02a3dfa9ae3a4969b29d4c0ecfaa6ad9
