#ifndef LF_USCALFE_FE_SPACE_LAGRANGE_O2_H
#define LF_USCALFE_FE_SPACE_LAGRANGE_O2_H

#include "uniform_scalar_fe_space.h"

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Quadratic Lagrangian Finite Element space
 *
 * Just a specialization of UniformScalarFESpace based on
 * FeLagrangeO2Tria, FeLagrangeO2Quad and FeLagrangePoint.
 */
template <typename SCALAR>
class FeSpaceLagrangeO2 : public UniformScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  /** @brief no default constructors */
  FeSpaceLagrangeO2() = delete;
  FeSpaceLagrangeO2(const FeSpaceLagrangeO2 &) = delete;
  FeSpaceLagrangeO2(FeSpaceLagrangeO2 &&) noexcept = default;
  FeSpaceLagrangeO2 &operator=(const FeSpaceLagrangeO2 &) = delete;
  FeSpaceLagrangeO2 &operator=(FeSpaceLagrangeO2 &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   */
  explicit FeSpaceLagrangeO2(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p)
      : UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<FeLagrangeO2Tria<SCALAR>>(),
            std::make_shared<FeLagrangeO2Quad<SCALAR>>(),
            std::make_shared<FeLagrangeO2Segment<SCALAR>>(),
            std::make_shared<FeLagrangePoint<SCALAR>>(2)) {}
  ~FeSpaceLagrangeO2() override = default;
};

}  // namespace lf::uscalfe

#endif  // LF_USCALFE_FE_SPACE_LAGRANGE_O2_H
