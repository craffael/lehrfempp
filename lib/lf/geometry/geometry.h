#ifndef __02a3dfa9ae3a4969b29d4c0ecfaa6ad9
#define __02a3dfa9ae3a4969b29d4c0ecfaa6ad9

/**
 * \brief Defines the Geometry interface and provides a number of
 *        classes that implement this interface + additional geometry related
 *        helper routines
 */
namespace lf::geometry {}

#include "geometry_interface.h"
#include "point.h"
#include "quad_o1.h"
#include "segment_o1.h"
#include "tria_o1.h"

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
enum struct RefinementPattern: int {
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

#endif  // __02a3dfa9ae3a4969b29d4c0ecfaa6ad9
