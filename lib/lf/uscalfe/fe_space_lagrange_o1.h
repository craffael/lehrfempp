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

#include "scalar_uniform_fe_space.h"

namespace lf::uscalfe {
/**
 * @brief Linear Lagrangian Finite Element space
 *
 * Just a specialization of ScalarUniformFESpace based on
 * FeLagrangeO1Tria, FeLagrangeO1Quad.
 *
 */
template <typename SCALAR>
class FeSpaceLagrangeO1 : public ScalarUniformFESpace<SCALAR> {
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
      : ScalarUniformFESpace<SCALAR>(
            mesh_p, std::make_shared<FeLagrangeO1Tria<SCALAR>>(),
            std::make_shared<FeLagrangeO1Quad<SCALAR>>(),
            std::make_shared<FeLagrangeO1Segment<SCALAR>>()) {}
  virtual ~FeSpaceLagrangeO1() = default;
};
}  // namespace lf::uscalfe

#endif  // __7a41e223dd0e4176af0371c2b57d2b67
