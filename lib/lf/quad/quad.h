/**
 * @file
 * @brief Main include file for the lf.quad module.
 * @author Raffael Casagrande
 * @date   2018-08-25 04:05:34
 * @copyright MIT License
 */

#ifndef INCG1f1e78e74a804e5490813ec6a9148231
#define INCG1f1e78e74a804e5490813ec6a9148231

#include "make_quad_rule.h"
#include "make_quad_rule_nodal.h"
#include "quad_rule_cache.h"

/** @brief Rules for numerical quadrature on reference entity shapes
 *
 * Quadrature rules are specified by
 * - a sequence of _reference coordinates_ for quadrature points
 * - a corresponding sequence of real-valued quadrature weights
 * The main interface is provided by the class @ref lf::quad::QuadRule.
 *
 * Refer to [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{par:lfquad}.
 *
 */
namespace lf::quad {}

#endif  // INCG1f1e78e74a804e5490813ec6a9148231
