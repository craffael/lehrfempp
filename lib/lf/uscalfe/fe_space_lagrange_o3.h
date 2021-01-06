/**
 * @file
 * @brief Defines a special FESpaceUniform that is optimized for third order
 *        lagrangian shape functions.
 * @author Philippe Peter
 * @date   November 2019
 * @copyright MIT License
 */

#ifndef LF_USCALFE_FE_SPACE_LAGRANGE_O3_H
#define LF_USCALFE_FE_SPACE_LAGRANGE_O3_H

#include "lagr_fe.h"
#include "uniform_scalar_fe_space.h"

namespace lf::uscalfe {
/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Cubic Lagrangian Finite Element space
 *
 * Just a specialization of UniformScalarFESpace based on
 * FeLagrangeO3Tria, FeLagrangeO3Quad, FeLagrangeO3Segment and FeLagrangePoint.
 *
 */
template <typename SCALAR>
class FeSpaceLagrangeO3 : public UniformScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  /** @brief no default constructors */
  FeSpaceLagrangeO3() = delete;
  FeSpaceLagrangeO3(const FeSpaceLagrangeO3 &) = delete;
  FeSpaceLagrangeO3(FeSpaceLagrangeO3 &&) noexcept = default;
  FeSpaceLagrangeO3 &operator=(const FeSpaceLagrangeO3 &) = delete;
  FeSpaceLagrangeO3 &operator=(FeSpaceLagrangeO3 &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   */
  explicit FeSpaceLagrangeO3(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p)
      : UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<FeLagrangeO3Tria<SCALAR>>(),
            std::make_shared<FeLagrangeO3Quad<SCALAR>>(),
            std::make_shared<FeLagrangeO3Segment<SCALAR>>(),
            std::make_shared<fe::FePoint<SCALAR>>(3)) {}
  ~FeSpaceLagrangeO3() override = default;
};
}  // namespace lf::uscalfe

#endif  // LF_USCALFE_FE_SPACE_LAGRANGE_O3_H
