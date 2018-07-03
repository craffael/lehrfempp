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
 * A refinement pattern has to be passed as first argument to the member function
 * ChildGeometry of a geometry object.
 *
 */
enum struct RefinementPattern: int {
  rp_copy, /**< just copy the geometry of the current entity */
  rp_split, /**< for edges: split into two halves */
  rp_bisec_0, /**< single bisection, refinement edge 0 */
  rp_bisec_1, /**< single bisection, refinement edge 1 */
  rp_bisec_2, /**< single bisection, refinement edge 2 */
  rp_trisect_01, /**< double bisection, splitting edges 0 and 1 */
  rp_trisect_02, /**< double bisection, splitting edges 0 and 2 */
  rp_trisect_10, /**< double bisection, splitting edges 1 and 0 */
  rp_trisect_12, /**< double bisection, splitting edges 1 and 2 */
  rp_trisect_20, /**< double bisection, splitting edges 2 and 0 */
  rp_trisect_21, /**< double bisection, splitting edges 2 and 1 */
  rp_regular, /**< standard "full" (regular) refinement */
};

#endif  // __02a3dfa9ae3a4969b29d4c0ecfaa6ad9
