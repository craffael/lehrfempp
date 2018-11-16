#ifndef LF_FESPACE_H
#define LF_FESPACE_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structure describing finite element spaces
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include "lagr_fe.h"

namespace lf::fe {

/**
 * @brief Space of scalar valued finite element functions on a hybrid 2D mesh
 *
 * The abstract concept of a (parametric) finite element space involves
 * - an underlying mesh
 * - the definition of a local set of shape functions for every cell.
 *   This class is restricted to parametric finite element spaces featuring
 *   the _same_ set of reference shape functions for every cell of a particular
 *   topological type. This is indicated by the attribute *Uniform* in the class
 *   name, cf., the class lf::assemble::UniformFEDofHandler.
 *
 * This class just contains (pointers to) objects representing the various building 
 * blocks of a finite element space. It does not offer elaborate methods. 
 */
class UniformScalarFiniteElementSpace {
 public:
  /** @brief default constructors, needed by std::vector
   * @note creates an invalid object that cannot be used. */
  UniformScalarFiniteElementSpace() = default;
  UniformScalarFiniteElementSpace(const UniformScalarFiniteElementSpace &) =
      delete;
  UniformScalarFiniteElementSpace(UniformScalarFiniteElementSpace &&) noexcept =
      default;
  UniformScalarFiniteElementSpace &operator=(
      const UniformScalarFiniteElementSpace &) = delete;
  UniformScalarFiniteElementSpace &operator=(
      UniformScalarFiniteElementSpace &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   * @param rfs_tria_p pointer to layout description for reference shape
   * functions on triangular cells
   * @param rfs_quad_p pointer to layout description for reference shape
   * functions on quadrilateral cells
   *
   * The schemes for local shape have to satisfy certain _compatibility
   * conditions_:
   * - nodes may carry at most one local/global shape function
   * - The number of interior shape functions for edges of triangles and
   * quadrilaterals must agree
   */
  UniformScalarFiniteElementSpace(
      std::shared_ptr<const lf::mesh::Mesh> mesh_p,
      std::unique_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
      std::unique_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p);

  /** @brief acess to underlying mesh
   *  @return a shared _pointer_ to the mesh
   */
  std::shared_ptr<const lf::mesh::Mesh> Mesh() const {
    LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
    return mesh_p_;
  }
  /** @brief acess to assciated local-to-global map
   * @return a reference to the lf::assemble::DofHandler object (immutable)
   */
  const lf::assemble::DofHandler &LocGlobMap() const {
    LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
    LF_VERIFY_MSG(dofh_p_ != nullptr,
                  "No valid FE space object: no dof handler");
    return *dofh_p_;
  }
  /** @brief access to layout of reference shape functions
   * @return a pair of references to ScalarReferenceFiniteElement<double>
   * objects. The first gives the reference shape function layou for triangular
   * cells, the second for quadrilateral cells.
   */
  std::pair<const ScalarReferenceFiniteElement<double> &,
            const ScalarReferenceFiniteElement<double> &>
  ShapeFunctionLayouts() const;
  
  /** @brief access to shape function layout for triangular cells */
  const ScalarReferenceFiniteElement<double> &TriaShapeFunctionLayout() const;
  
  /** @brief access to shape function layout for quadrilateral cells */
  const ScalarReferenceFiniteElement<double> &QuadShapeFunctionLayout() const;

  /** @brief number of _interior_ shape functions associated to entities of various types
   */
  size_type NumRefShapeFunctions(lf::base::RefEl ref_el_type) const;
  
  /** @brief No special destructor */
  virtual ~UniformScalarFiniteElementSpace() = default;
  
 private:
  /** Underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  /** Description of reference shape functions on triangular cells */
  std::unique_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p_;
  /** Description of reference shape functions on quadrilateral cells */
  std::unique_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p_;
  /** Local-to-global index map for the finite element space */
  std::unique_ptr<lf::assemble::UniformFEDofHandler> dofh_p_;
  /** Checks whether some pointer are not valid */
  bool check_ptr() const {
    LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
    LF_VERIFY_MSG(dofh_p_ != nullptr,
                  "No valid FE space object: no dof handler");
    LF_VERIFY_MSG((rfs_quad_p_ != nullptr) && (rfs_quad_p_ != nullptr),
                  "No valid FE space object: no rsfs");
    return true;
  }
 public:
  /** Output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_mesh = 1;
  static const unsigned int kout_dofh = 2;
  static const unsigned int kout_rsfs = 4;
};  // end class definition UniformScalarFiniteElementSpace

  /** @brief output operator for scalar parametric finite element space */
  std::ostream &operator<<(std::ostream &o, const UniformScalarFiniteElementSpace &fes);
  
}  // namespace lf::fe

#endif
