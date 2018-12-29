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
 * This class just contains (pointers to) objects representing the various
 * building blocks of a finite element space. It does not offer elaborate
 * methods.
 *
 * @note Some of the pointers may be NULL. For instance, if all computations
 *       are done on purely triangular meshes then a finite element
 * specification for quadrilaterals need not be given, \see
 * FeSpaceUniformScalar().
 */
template <typename SCALAR>
class FeSpaceUniformScalar {
 public:
  /** @brief default constructors, needed by std::vector
   * @note creates an invalid object that cannot be used. */
  FeSpaceUniformScalar() = delete;
  FeSpaceUniformScalar(const FeSpaceUniformScalar &) = delete;
  FeSpaceUniformScalar(FeSpaceUniformScalar &&) noexcept = default;
  FeSpaceUniformScalar &operator=(const FeSpaceUniformScalar &) = delete;
  FeSpaceUniformScalar &operator=(FeSpaceUniformScalar &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   * @param rfs_tria_p pointer to layout description for reference shape
   * functions on triangular cells
   * @param rfs_quad_p pointer to layout description for reference shape
   * functions on quadrilateral cells
   * @param rfs_edge_p pointer to layout description for reference shape
   * functions on the edges
   *
   * The schemes for local shape have to satisfy certain _compatibility
   * conditions_:
   * - nodes may carry at most one local/global shape function
   * - The number of interior shape functions for edges of triangles and
   * quadrilaterals must agree.
   *
   * @note none of the shape function layouts needs to be specified; just pass
   *       a null pointer. This will then restrict the applicability of
   *       the resulting finite element space objects to particular meshes.
   */
  FeSpaceUniformScalar(
      std::shared_ptr<const lf::mesh::Mesh> mesh_p,
      std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> rfs_tria_p,
      std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> rfs_quad_p,
      std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> rfs_edge_p =
          nullptr)
      : mesh_p_(mesh_p),
        rfs_tria_p_(rfs_tria_p),
        rfs_quad_p_(rfs_quad_p),
        rfs_edge_p_(rfs_edge_p) {
    init();
  }

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

  /** @brief access to shape function layout for cells
   *
   * @param ref_el_type type of entit, can be anything except for
   * lf::base::RefEl::kPoint()
   *
   * @note NULL pointers may be returned by this method in case a finite element
   * specification was not given for a particular topological type of entity.
   */
  std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>>
  ShapeFunctionLayout(lf::base::RefEl rel_el_type) const;

  /** @brief number of _interior_ shape functions associated to entities of
   * various types
   */
  size_type NumRefShapeFunctions(lf::base::RefEl ref_el_type) const;

  /** @brief No special destructor */
  virtual ~FeSpaceUniformScalar() = default;

 private:
  /** Underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  /** Description of reference shape functions on triangular cells */
  std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> rfs_tria_p_;
  /** Description of reference shape functions on quadrilateral cells */
  std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> rfs_quad_p_;
  /** Description of reference shape functions on an edge */
  std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> rfs_edge_p_;
  /** Numbers of local shape functions for different types of entities */
  size_type num_rsf_node_{0}, num_rsf_edge_{0}, num_rsf_tria_{0},
      num_rsf_quad_{0};
  /** Local-to-global index map for the finite element space */
  std::unique_ptr<lf::assemble::UniformFEDofHandler> dofh_p_;

  /** Initialization of class member variables and consistency checks */
  void init(void);
  /** Checks whether some pointer are not valid */
  bool check_ptr() const {
    LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
    LF_VERIFY_MSG(dofh_p_ != nullptr,
                  "No valid FE space object: no dof handler");
    LF_VERIFY_MSG((rfs_quad_p_ != nullptr) && (rfs_quad_p_ != nullptr),
                  "No valid FE space object: no rsfs for cells");
    return true;
  }

 public:
  /** Output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_mesh = 1;
  static const unsigned int kout_dofh = 2;
  static const unsigned int kout_rsfs = 4;
};  // end class definition FeSpaceUniformScalar

/** @brief output operator for scalar parametric finite element space */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o,
                         const FeSpaceUniformScalar<SCALAR> &fes);

// Output control variable
template <typename SCALAR>
unsigned int FeSpaceUniformScalar<SCALAR>::ctrl_ = 0;

// Initialization methods
template <typename SCALAR>
void FeSpaceUniformScalar<SCALAR>::init() {
  LF_VERIFY_MSG(mesh_p_ != nullptr, "Missing mesh!");
  LF_VERIFY_MSG((rfs_quad_p_ != nullptr) || (rfs_tria_p_ != nullptr),
                "Missing FE specification for cells");
  LF_VERIFY_MSG((mesh_p_->DimMesh() == 2), "Only for 2D meshes");

  // Check whether all required finite element specifications are provided
  LF_VERIFY_MSG((mesh_p_->NumEntities(lf::base::RefEl::kTria()) == 0) ||
                    (rfs_tria_p_ != nullptr),
                "Missing FE specification for triangles");
  LF_VERIFY_MSG((mesh_p_->NumEntities(lf::base::RefEl::kQuad()) == 0) ||
                    (rfs_quad_p_ != nullptr),
                "Missing FE specification for quads");

  // Compatibility checks and initialization of numbers of shape functions
  // In particular only a single shape function may be associated to a node
  if (rfs_tria_p_ != nullptr) {
    LF_VERIFY_MSG((*rfs_tria_p_).RefEl() == lf::base::RefEl::kTria(),
                  "Wrong type for triangle!");
    LF_VERIFY_MSG((*rfs_tria_p_).NumRefShapeFunctions(2) <= 1,
                  "At most one shape function can be assigned to each vertex");
    num_rsf_node_ = (*rfs_tria_p_).NumRefShapeFunctions(2);
    num_rsf_edge_ = (*rfs_tria_p_).NumRefShapeFunctions(1);
    num_rsf_tria_ = (*rfs_tria_p_).NumRefShapeFunctions(0);
  }
  if (rfs_quad_p_ != nullptr) {
    LF_VERIFY_MSG((*rfs_quad_p_).RefEl() == lf::base::RefEl::kQuad(),
                  "Wrong type for quad!");
    LF_VERIFY_MSG((*rfs_quad_p_).NumRefShapeFunctions(2) <= 1,
                  "At most one shape function can be assigned to each vertex");
    num_rsf_node_ = (*rfs_quad_p_).NumRefShapeFunctions(2);
    num_rsf_edge_ = (*rfs_quad_p_).NumRefShapeFunctions(1);
    num_rsf_quad_ = (*rfs_quad_p_).NumRefShapeFunctions(0);
  }
  if (rfs_edge_p_ != nullptr) {
    LF_VERIFY_MSG((*rfs_edge_p_).RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for edge!");
    LF_VERIFY_MSG((*rfs_edge_p_).NumRefShapeFunctions(1) <= 1,
                  "At most one shape function can be assigned to each vertex");
    num_rsf_node_ = (*rfs_edge_p_).NumRefShapeFunctions(1);
    num_rsf_edge_ = (*rfs_edge_p_).NumRefShapeFunctions(0);
  }

  // Compatibility check for numbers of local shape functions associated with
  // edges
  if ((rfs_tria_p_ != nullptr) && (rfs_quad_p_ != nullptr)) {
    LF_ASSERT_MSG(((*rfs_tria_p_).NumRefShapeFunctions(2) ==
                   (*rfs_quad_p_).NumRefShapeFunctions(2)),
                  "#RSF mismatch on nodes "
                      << (*rfs_tria_p_).NumRefShapeFunctions(2) << " <-> "
                      << (*rfs_quad_p_).NumRefShapeFunctions(2));
    LF_ASSERT_MSG(((*rfs_tria_p_).NumRefShapeFunctions(1) ==
                   (*rfs_quad_p_).NumRefShapeFunctions(1)),
                  "#RSF mismatch on edges "
                      << (*rfs_tria_p_).NumRefShapeFunctions(1) << " <-> "
                      << (*rfs_quad_p_).NumRefShapeFunctions(1));
  }
  if ((rfs_tria_p_ != nullptr) && (rfs_edge_p_ != nullptr)) {
    LF_ASSERT_MSG(((*rfs_tria_p_).NumRefShapeFunctions(2) ==
                   (*rfs_edge_p_).NumRefShapeFunctions(1)),
                  "#RSF mismatch on nodes "
                      << (*rfs_tria_p_).NumRefShapeFunctions(2) << " <-> "
                      << (*rfs_edge_p_).NumRefShapeFunctions(1));
    LF_ASSERT_MSG(((*rfs_tria_p_).NumRefShapeFunctions(1) ==
                   (*rfs_edge_p_).NumRefShapeFunctions(0)),
                  "#RSF mismatch on edges "
                      << (*rfs_tria_p_).NumRefShapeFunctions(1) << " <-> "
                      << (*rfs_edge_p_).NumRefShapeFunctions(0));
  }
  if ((rfs_quad_p_ != nullptr) && (rfs_edge_p_ != nullptr)) {
    LF_ASSERT_MSG(((*rfs_quad_p_).NumRefShapeFunctions(2) ==
                   (*rfs_edge_p_).NumRefShapeFunctions(1)),
                  "#RSF mismatch on edges "
                      << (*rfs_quad_p_).NumRefShapeFunctions(2) << " <-> "
                      << (*rfs_edge_p_).NumRefShapeFunctions(1));
    LF_ASSERT_MSG(((*rfs_quad_p_).NumRefShapeFunctions(1) ==
                   (*rfs_edge_p_).NumRefShapeFunctions(0)),
                  "#RSF mismatch on edges "
                      << (*rfs_quad_p_).NumRefShapeFunctions(1) << " <-> "
                      << (*rfs_edge_p_).NumRefShapeFunctions(0));
  }

  // Initialization of dof handler starting with collecting the number of
  // interior reference shape functions
  lf::assemble::UniformFEDofHandler::dof_map_t rsf_layout{
      {lf::base::RefEl::kPoint(), num_rsf_node_},
      {lf::base::RefEl::kSegment(), num_rsf_edge_},
      {lf::base::RefEl::kTria(), num_rsf_tria_},
      {lf::base::RefEl::kQuad(), num_rsf_quad_}};
  dofh_p_ =
      std::make_unique<lf::assemble::UniformFEDofHandler>(mesh_p_, rsf_layout);
}

template <typename SCALAR>
std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>>
FeSpaceUniformScalar<SCALAR>::ShapeFunctionLayout(
    lf::base::RefEl ref_el_type) const {
  // Retrieve specification of local shape functions
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      LF_VERIFY_MSG(false, "No shape functions for point!");
      break;
    }
    case lf::base::RefEl::kSegment(): {
      // LF_ASSERT_MSG(rfs_edge_p_ != nullptr, "No RSF for edges!");
      return rfs_edge_p_;
    }
    case lf::base::RefEl::kTria(): {
      // LF_ASSERT_MSG(rfs_tria_p_ != nullptr, "No RSF for triangles!");
      return rfs_tria_p_;
    }
    case lf::base::RefEl::kQuad(): {
      // LF_ASSERT_MSG(rfs_quad_p_ != nullptr, "No RSF for quads!");
      return rfs_quad_p_;
    }
    default: { LF_VERIFY_MSG(false, "Illegal entity type"); }
  }
  return nullptr;
}

/* number of _interior_ shape functions associated to entities of various types
 */
template <typename SCALAR>
size_type FeSpaceUniformScalar<SCALAR>::NumRefShapeFunctions(
    lf::base::RefEl ref_el_type) const {
  LF_ASSERT_MSG((rfs_quad_p_ != nullptr) && (rfs_quad_p_ != nullptr),
                "No valid FE space object: no rsfs");
  // Retrieve number of interior shape functions from rsf layouts
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return num_rsf_node_;
    }
    case lf::base::RefEl::kSegment(): {
      return num_rsf_edge_;
    }
    case lf::base::RefEl::kTria(): {
      return num_rsf_tria_;
    }
    case lf::base::RefEl::kQuad(): {
      return num_rsf_quad_;
    }
    dafault : { LF_VERIFY_MSG(false, "Illegal entity type"); }
  }
  return 0;
}

/** output operator for scalar parametric finite element space */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o,
                         const FeSpaceUniformScalar<SCALAR> &fes) {
  o << "Uniform scalar FE space, dim = " << fes.LocGlobMap().NoDofs()
    << std::endl;
  if (FeSpaceUniformScalar<SCALAR>::ctrl_ &
      FeSpaceUniformScalar<SCALAR>::kout_mesh) {
    o << fes.Mesh() << std::endl;
  }
  if (FeSpaceUniformScalar<SCALAR>::ctrl_ &
      FeSpaceUniformScalar<SCALAR>::kout_dofh) {
    o << fes.LocGlobMap() << std::endl;
  }
  if (FeSpaceUniformScalar<SCALAR>::ctrl_ &
      FeSpaceUniformScalar<SCALAR>::kout_rsfs) {
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

}  // namespace lf::fe

#endif
