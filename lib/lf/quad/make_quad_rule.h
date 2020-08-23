/**
 * @file
 * @brief Declaration of various routines for creating QuadRule objects
 * @author Raffael Casagrande
 * @date   2018-08-19 06:29:37
 * @copyright MIT License
 */

#include <lf/base/base.h>
#include "quad_rule.h"

namespace lf::quad {

/**
 * @brief Returns a QuadRule object for the given Reference Element and Degree
 * @param ref_el The type of reference element
 * @param degree The minimum degree that the QuadRule object should have.
 * @return A QuadRule object for the given reference element and with a degree
 * >= `degree`
 *
 * This method tries to return optimal quadrature rules when possible:
 * - For Segments it returns Gauss-Legendre quadrature rules
 * - For Triangles, it returns specific, hard-coded quadrature rules up to order
 * 50, afterwards, the Duffy-Transform to map a tensor Gaussian Quadrature on a
 * Square to the triangle.
 * - For Quadrilaterals it uses tensor products of Gauss-Legendre rules
 */
QuadRule make_QuadRule(base::RefEl ref_el, unsigned degree);

/** @defgroup namedqr Special "named" quadrature rules
 * @breif Creation of special quadrature rules
 *
 * These functions provide a number of special quadrature rules
 * in the forms of @ref QuadRule objects. They are meant to be
 * used in low-order finite element methods when the quadrature
 * rules can easily fixed a-priori and should feature a special
 * location of quadrature points as is required, e.g., for mass-lumping
 * techniques.
 *
 * For some examples see [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{ex:triquadrules}.
 *
 * QuadRule objects with a particular order but otherwise unspecified properties
 * can be created by @ref lf::quad::make_QuadRule().
 * @{
 */

/** @brief midpoint quadrature rule for triangles
 *
 * This quadrature rule relies on the center of gravity as single quadrature
 * node
 */
QuadRule make_TriaQR_MidpointRule();
inline QuadRule make_TriaQR_P1O2() { return make_TriaQR_MidpointRule(); }

/** @brief edge midpoint quadrature rule for reference triangles
 *
 * This quadrature rule relies on point evaluations at the midpoints of
 * all edges with equal weights = 1/6 for the reference triangle.
 *
 * The rule is exact for quadratic bi-variate polynomials
 */
QuadRule make_TriaQR_EdgeMidpointRule();
inline QuadRule make_TriaQR_P3O3() { return make_TriaQR_EdgeMidpointRule(); }

/** @brief Seven point triangular quadrature rule of order 6 */
QuadRule make_TriaQR_P7O6();

/** @brief Six point triangular quadrature rule of order 4 */
QuadRule make_TriaQR_P6O4();

/** @brief midpoint quadrature rule for quadrilaterals
 *
 * This quadrature rule relies on the center of gravity as single quadrature
 * node. It agrees with the lowest-order tensor product Gauss rule
 */
inline QuadRule make_QuadQR_MidpointRule() {
  return make_QuadRule(lf::base::RefEl::kQuad(), 1);
}
inline QuadRule make_QuadQR_P1O2() { return make_QuadQR_MidpointRule(); }

/** @brief edge midpoint quadrature rule for unit square (= reference quad)
 *
 * This quadrature rule relies on point evaluations at the midpoints of
 * all edges with equal weights = 1/4 for the unit square.
 *
 * The rule is exact for bilinear polynomials.
 */
QuadRule make_QuadQR_EdgeMidpointRule();
inline QuadRule make_QuadQR_P4O2() { return make_QuadQR_EdgeMidpointRule(); }

/** @brief Fourth-order tensor product Gauss rule for quadrilaterals
 *
 * Rule uses four interior quadrature nodes and is exact for tensor-product
 * polynomials up to degree 3.
 */
inline QuadRule make_QuadQR_P4O4() {
  return make_QuadRule(lf::base::RefEl::kQuad(), 3);
}

/** @} */
}  // namespace lf::quad
