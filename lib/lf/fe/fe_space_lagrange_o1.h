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

#include "fe_space_uniform_scalar.h"

namespace lf::fe {
/**
 * @brief Linear Lagrangian Finite Element space
 *
 * Just a specialization of FeSpaceUniformScalar based on
 * TriaLinearLagrangeFE, QuadLinearLagrangeFE.
 *
 */
template <typename SCALAR>
class FeSpaceLagrangeO1 : public FeSpaceUniformScalar<SCALAR> {
 public:
  /** @brief default constructors, needed by std::vector
   * @note creates an invalid object that cannot be used. */
  FeSpaceLagrangeO1() = default;
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
      : FeSpaceUniformScalar<SCALAR>(
            mesh_p, std::make_shared<TriaLinearLagrangeFE<SCALAR>>(),
            std::make_shared<QuadLinearLagrangeFE<SCALAR>>(),
            std::make_shared<SegmentLinearLagrangeFE<SCALAR>>()) {}
  virtual ~FeSpaceLagrangeO1() {}
};
}  // namespace lf::fe

#endif  // __7a41e223dd0e4176af0371c2b57d2b67
