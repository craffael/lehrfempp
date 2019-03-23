/**
 * @file
 * @brief Declaration of ChildGeometry() tests for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-30 18:08:17
 * @copyright MIT License
 */

#ifndef LEHRFEMPP_CHECK_CHILD_GEOMETRY_H
#define LEHRFEMPP_CHECK_CHILD_GEOMETRY_H

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

namespace lf::geometry::test_utils {

/**
 * @brief Checks that the mapping geom.ChildGeometry() is the same as geom
 *    composed with the mapping imposed by ref_pat
 * @param geom The geometry object whose ChildGeometry() method should be
 * checked
 * @param ref_pat The refinement pattern that is used for the test
 * @param qr_provider Provides the quadrature rules whose points are used to
 *   check whether the two mappings agree
 *
 * @warning For non-linear mappings, this test can fail, especially for
 *   quadrilaterals which are split into triangles! In this case it makes sense
 *   to use a nodal quadrature rule for the test (see
 *   quad::make_QuadRuleNodal())
 */
void checkChildGeometry(
    const lf::geometry::Geometry &geom,
    const lf::geometry::RefinementPattern &ref_pat,
    const std::function<lf::quad::QuadRule(lf::base::RefEl)> &qr_provider);

/**
 * @brief Checks if the total volume is conserved after call to
 * ChildGeometry()
 * @param geom geometry object to be evaluated
 * @param refPat refinement pattern for child generation
 * @param anchor local number of anchor element for refinement pattern
 */
void checkChildGeometryVolume(
    const lf::geometry::Geometry &geom, const lf::refinement::RefPat &refPat,
    const lf::base::sub_idx_t &anchor = lf::refinement::idx_nil);

}  // namespace lf::geometry::test_utils

#endif  // LEHRFEMPP_CHECK_CHILD_GEOMETRY_H
