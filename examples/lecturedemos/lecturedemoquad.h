#ifndef LF_LD_QUAD_H
#define LF_LD_QUAD_H
/**
 * @file
 * @brief Simple LehrFEM++ demo and sample codes for the use of numerical
 * quadrature
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/quad/quad.h>
#include "lecturedemo.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

namespace lecturedemo {
/**
 * @brief Entity-oriented numerical integration of a function
 *
 * @tparm FUNCTOR functor type for scalar-valued function
 * @param mesh reference to a mesh object
 * @param quadrules associative array of quadrature rules, key is the entity
 * type
 * @param f functor object providing the function to be integrated
 * @param codim co-dimension of entities to be visited
 * @param pred selector for entities that should be covered by integration
 *
 * This functions performs entity-based composite numerical quadrature of a
 * function \f$f:\Omega\mapsto V\f$, where \f$V\f$ is a vector space. The domain
 * of integration is the union of all mesh entities of co-dimension `codim` for
 * which the predicate `pred` returns `true`. The functin uses the local
 * quadrature rule(s) passed through `quadrules`.
 *
 * ##Requirements for FUNCTOR
 *
 * A FUNCTOR type object must provide `V operator () (Eigen::VectorXd x)`, where
 * `V` is a type that fits the concept of a real vector space.
 */
//clang-format off
/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
auto localQuadFunction(const lf::mesh::Mesh &mesh,
                       std::map<lf::base::RefEl, lf::quad::QuadRule> quadrules,
                       FUNCTOR &&f, dim_t codim,
                       std::function<bool(const lf::mesh::Entity &)> pred =
                           [](const lf::mesh::Entity & /*entity*/) -> bool {
                         return true;
                       }) {
  LF_ASSERT_MSG(mesh.DimMesh() >= codim, "Illegal codim = " << codim);
  // Variable for summing the result
  using value_t = std::invoke_result_t<FUNCTOR, Eigen::VectorXd>;
  value_t sum_var{};
  // Loop over entities of co-dimension codim
  for (const lf::mesh::Entity &entity : mesh.Entities(codim)) {
    // Obtain geometry information for entity
    const lf::geometry::Geometry &geo{*entity.Geometry()};
    // obtain quadrature rule suitable for entity type
    auto tmp = quadrules.find(entity.RefEl());
    if (tmp != quadrules.end()) {
      // A quadrature rule has been found
      const lf::quad::QuadRule &qr{tmp->second};
      // Number of quadrature points
      const size_type P = qr.NumPoints();
      // Quadrature points
      const Eigen::MatrixXd zeta_ref{qr.Points()};
      // Map quadrature points to physical/world coordinates
      const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
      // Quadrature weights
      const Eigen::VectorXd w_ref{qr.Weights()};
      // Gramian determinants
      const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};
      // Iterate over the quadrature points
      for (int l = 0; l < P; ++l) {
        sum_var += w_ref[l] * f(zeta.col(l)) * gram_dets[l];
      }
    } else {
      LF_VERIFY_MSG(false, "Missing quadrature rule for " << entity.RefEl());
    }
  }
  return sum_var;
}

/* SAM_LISTING_END_1 */
// clang-format on

/** @brief Driver function for quadrature demo */
void lecturedemoquad();

}  // namespace lecturedemo

#endif
