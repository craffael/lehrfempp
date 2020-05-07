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

#include "lagr_fe.h"

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
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
 * This class is covered in @\lref{par:fespace}.
 */
template <typename SCALAR>
class ScalarFESpace {
 public:
  using Scalar = SCALAR;

  /** @brief default constructor, needed by std::vector
   * @note creates an invalid object that cannot be used. */
  ScalarFESpace() = default;
  ScalarFESpace(const ScalarFESpace &) = delete;
  ScalarFESpace(ScalarFESpace &&) noexcept = default;
  ScalarFESpace &operator=(const ScalarFESpace &) = delete;
  ScalarFESpace &operator=(ScalarFESpace &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   */
  ScalarFESpace(std::shared_ptr<const lf::mesh::Mesh> mesh_p)
      : mesh_p_(std::move(mesh_p)) { }

  /** @brief acess to underlying mesh
   *  @return a shared _pointer_ to the mesh
   */
  [[nodiscard]] std::shared_ptr<const lf::mesh::Mesh> Mesh() const {
    LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
    return mesh_p_;
  }
  /** @brief access to associated local-to-global map
   * @return a reference to the lf::assemble::DofHandler object (immutable)
   */
  [[nodiscard]] virtual const lf::assemble::DofHandler &LocGlobMap() const = 0;

  /** @brief access to shape function layout for cells
   *
   * @param entity The entity to get the reference element for
   *
   * @warning NULL pointers may be returned by this method in case a finite
   * element specification was not given for a particular topological type of
   * entity.
   */
  [[nodiscard]] virtual std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>>
  ShapeFunctionLayout(const lf::mesh::Entity &entity) const = 0;

  /** @brief number of _interior_ shape functions associated to entities of
   * various types
   */
  [[nodiscard]] virtual size_type NumRefShapeFunctions(const lf::mesh::Entity &entity) const = 0;

  /** @brief No special destructor */
  virtual ~ScalarFESpace() = default;

 private:
  /** Underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;

 public:
  /** Output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_mesh = 1;
  static const unsigned int kout_dofh = 2;
  static const unsigned int kout_rsfs = 4;
};  // end class definition ScalarFESpace

/** @brief output operator for scalar parametric finite element space */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o,
                         const ScalarFESpace<SCALAR> &fes);

// Output control variable
template <typename SCALAR>
unsigned int ScalarFESpace<SCALAR>::ctrl_ = 0;

/** output operator for scalar parametric finite element space */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o,
                         const ScalarFESpace<SCALAR> &fes) {
  o << "Uniform scalar FE space, dim = " << fes.LocGlobMap().NumDofs()
    << std::endl;
  if (ScalarFESpace<SCALAR>::ctrl_ &
      ScalarFESpace<SCALAR>::kout_mesh) {
    o << fes.Mesh() << std::endl;
  }
  if (ScalarFESpace<SCALAR>::ctrl_ &
      ScalarFESpace<SCALAR>::kout_dofh) {
    o << fes.LocGlobMap() << std::endl;
  }
  if (ScalarFESpace<SCALAR>::ctrl_ &
      ScalarFESpace<SCALAR>::kout_rsfs) {
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kPoint()) << " rsfs @ nodes"
      << std::endl;
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kSegment())
      << " rsfs @ edges" << std::endl;
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kTria())
      << " rsfs @ triangles" << std::endl;
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kQuad()) << " rsfs @ quads"
      << std::endl;
  }
  return o;
}

}  // namespace lf::uscalfe

#endif
