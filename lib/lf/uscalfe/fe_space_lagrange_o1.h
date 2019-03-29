/**
 * @file
 * @brief Defines a special FESpaceUniform that is optimized for first order
 *        lagrangian shape functions.
 * @author Raffael Casagrande
 * @date   2018-12-29 03:05:42
 * @copyright MIT License
 */

#ifndef __7a41e223dd0e4176af0371c2b57d2b67
#define __7a41e223dd0e4176af0371c2b57d2b67

#include "uniform_scalar_fe_space.h"

namespace lf::uscalfe {
/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief (Bi)Linear Lagrangian Finite Element space
 *
 * Just a specialization of UniformScalarFESpace based on
 * FeLagrangeO1Tria, FeLagrangeO1Quad and FeLagrangePoint.
 *
 * ### Example
 * @snippet fe_space_lagrange_o1.cc usage
 *
 */
template <typename SCALAR>
class FeSpaceLagrangeO1 : public UniformScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  /** @brief no default constructors */
  FeSpaceLagrangeO1() = delete;
  FeSpaceLagrangeO1(const FeSpaceLagrangeO1 &) = delete;
  FeSpaceLagrangeO1(FeSpaceLagrangeO1 &&) noexcept = default;
  FeSpaceLagrangeO1 &operator=(const FeSpaceLagrangeO1 &) = delete;
  FeSpaceLagrangeO1 &operator=(FeSpaceLagrangeO1 &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   */
  explicit FeSpaceLagrangeO1(std::shared_ptr<const lf::mesh::Mesh> mesh_p)
      : UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<FeLagrangeO1Tria<SCALAR>>(),
            std::make_shared<FeLagrangeO1Quad<SCALAR>>(),
            std::make_shared<FeLagrangeO1Segment<SCALAR>>(),
            std::make_shared<FeLagrangePoint<SCALAR>>(1)) {}
  ~FeSpaceLagrangeO1() override = default;
};
}  // namespace lf::uscalfe

#endif  // __7a41e223dd0e4176af0371c2b57d2b67
