#ifndef LF_USCALFE_SCALAR_FE_SPACE_H_
#define LF_USCALFE_SCALAR_FE_SPACE_H_
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structure describing scalar-valued finite element spaces
 * @author Ralf Hiptmair, Tobias Rohner
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>

#include "scalar_reference_finite_element.h"

namespace lf::fe {

/**
 * @headerfile lf/fe/fe.h
 * @brief Space of scalar valued finite element functions on a \ref mesh::Mesh
 * "Mesh"
 *
 * @tparam SCALAR Scalar type of the finite element functions, usually either
 * `double` or `complex<double>`
 *
 * A ScalarFESpace can be thought of as a set of basis functions \f$ b^1, \ldots
 * b^N \f$ which span a (approximation) space \f$ V_h := \operatorname{span} \{
 * b^1, \ldots b^N\} \f$. The basis functions \f$ b^i\f$ are usually selected
 * such that their support is limited to a small number of mesh cells. This is
 * reflected in the design of ScalarFESpace:
 * - ShapeFunctionLayout() returns a \ref ScalarReferenceFiniteElement for an
 *   \ref mesh::Entity "entity" which describes only the shape functions which
 *   are non-zero on this entity (respectively whose support intersects with the
 *   entity). In other words, the ScalarReferenceFiniteElement for an entity \f$
 *   e \f$ describes the restriction of the approximation space
 *   \f$ \left. V_h \right|_e \f$ (keep in mind that for technical reasons the
 *   shape functions of the ScalarReferenceFiniteElement are defined on a \ref
 *   base::RefEl "reference element"!)
 * - LocGlobMap() returns a \ref assemble::DofHandler "DofHandler" which
 *   identifies the local basis function (returned by
 *   ShapeFunctionLayout()) with some of the global basis functions \f$
 *   b^1\ldots b^N \f$
 *
 * This class is covered in @lref{par:fespace}.
 *
 * #### Relation to uscalfe::UniformScalarFESpace
 * As you can see in the inheritance diagram, every
 * \ref uscalfe::UniformScalarFESpace "UniformScalarFESpace" is also a
 * ScalarFESpace. The UniformScalarFESpace class and all other classes in the
 * namespace `lf::uscalfe` make additionally the assumption that when we
 * restrict the basis functions \f$ b^1, \ldots b^N\f$ to a cell, the remaining
 * set of basis functions are always the same for every cell/reference element
 * (when defined on the \ref base::RefEl "reference element").
 * I.e. ShapeFunctionLayout() will always return the same
 * ScalarReferenceFiniteElement for the same \ref base::RefEl
 * "reference element".
 *
 */
template <typename SCALAR>
class ScalarFESpace {
 public:
  using Scalar = SCALAR;

 protected:
  /** @brief default constructor, needed by std::vector
   * @note creates an invalid object that cannot be used. */
  ScalarFESpace() = default;
  ScalarFESpace(const ScalarFESpace &) = default;
  ScalarFESpace(ScalarFESpace &&) noexcept = default;
  ScalarFESpace &operator=(const ScalarFESpace &) = default;
  ScalarFESpace &operator=(ScalarFESpace &&) noexcept = default;

 public:
  /** @brief acess to underlying mesh
   *  @return a shared _pointer_ to the mesh
   */
  [[nodiscard]] virtual std::shared_ptr<const lf::mesh::Mesh> Mesh() const = 0;
  /** @brief access to associated local-to-global map
   *
   * @return a reference to the lf::assemble::DofHandler object (immutable)
   */
  [[nodiscard]] virtual const lf::assemble::DofHandler &LocGlobMap() const = 0;

  /** @brief access to shape function layout for mesh entities
   *
   * @param entity The entity to get the reference element for
   *
   * @warning NULL pointers may be returned by this method in case a finite
   * element specification was not given for a particular topological type of
   * entity.
   *
   * @note The returned ShapeFunctionLayout pointer will remain for the entire
   * lifetime of the owning ScalarFESpace.
   *
   * @see ScalarReferenceFiniteElement
   */
  [[nodiscard]] virtual ScalarReferenceFiniteElement<SCALAR> const *
  ShapeFunctionLayout(const lf::mesh::Entity &entity) const = 0;

  /** @brief number of _interior_ shape functions associated to a particular
   * mesh entity.
   *
   * @param entity mesh entity to be queried
   * @return number of _interior_ shape functions
   */
  [[nodiscard]] virtual size_type NumRefShapeFunctions(
      const lf::mesh::Entity &entity) const = 0;

  /** @brief No special destructor */
  virtual ~ScalarFESpace() = default;
};  // end class definition ScalarFESpace

/** @brief output operator for scalar parametric finite element space */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o, const ScalarFESpace<SCALAR> &fes) {
  fes.PrintInfo(o, 0);
  return o;
}

}  // namespace lf::fe

#endif
