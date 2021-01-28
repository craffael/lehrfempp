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
 * @brief Space of scalar valued finite element functions on a _hybrid 2D mesh_
 *
 * @tparam SCALAR underlying scalar type, usually either `double` or
 * `complex<double>`
 *
 * The abstract concept of a (parametric) finite element space involves
 * - an underlying mesh
 * - the definition of a local set of shape functions for every cell.
 *
 * This class just contains (pointers to) objects representing the various
 * building blocks of a finite element space. It does not offer elaborate
 * methods.
 *
 * @note Some of the pointers may be NULL. For instance, if all computations
 *       are done on purely triangular meshes then a finite element
 * specification for quadrilaterals need not be given.
 *
 * This class is covered in @lref{par:fespace}.
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
