#ifndef LF_FESPACE_H
#define LF_FESPACE_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structure describing scalar-valued finite element spaces
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>

#include <lf/fe/scalar_fe_space.h>
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
 * specification for quadrilaterals need not be given.
 *
 * This class is covered in @\lref{par:fespace}.
 */
template <typename SCALAR>
class UniformScalarFESpace : public lf::fe::ScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  /** @brief default constructor, needed by std::vector
   * @note creates an invalid object that cannot be used. */
  UniformScalarFESpace() = default;
  UniformScalarFESpace(const UniformScalarFESpace &) = delete;
  UniformScalarFESpace(UniformScalarFESpace &&) noexcept = default;
  UniformScalarFESpace &operator=(const UniformScalarFESpace &) = delete;
  UniformScalarFESpace &operator=(UniformScalarFESpace &&) noexcept = default;
  /**
   * @brief Main constructor: sets up the local-to-global index mapping (dof
   * handler)
   *
   * @param mesh_p shared pointer to underlying mesh (immutable)
   * @param rfs_tria_p shared pointer to layout description for reference shape
   * functions on triangular cells
   * @param rfs_quad_p shared pointer to layout description for reference shape
   * functions on quadrilateral cells
   * @param rfs_edge_p shared pointer to layout description for reference shape
   * functions on the edges
   * @param rfs_point_p shared pointer to layout description for reference shape
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
  UniformScalarFESpace(
      std::shared_ptr<const lf::mesh::Mesh> mesh_p,
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
          rfs_tria_p,
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
          rfs_quad_p,
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
          rfs_edge_p = nullptr,
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
          rfs_point_p = nullptr)
      : lf::fe::ScalarFESpace<SCALAR>(),
        rfs_tria_p_(std::move(rfs_tria_p)),
        rfs_quad_p_(std::move(rfs_quad_p)),
        rfs_edge_p_(std::move(rfs_edge_p)),
        rfs_point_p_(std::move(rfs_point_p)) {
    init(std::move(mesh_p));
  }

  /** @brief acess to underlying mesh
   *  @return a shared _pointer_ to the mesh
   */
  [[nodiscard]] std::shared_ptr<const lf::mesh::Mesh> Mesh() const {
    return dofh_p_->Mesh();
  }

  /** @brief access to associated local-to-global map
   * @return a reference to the lf::assemble::DofHandler object (immutable)
   */
  [[nodiscard]] const lf::assemble::DofHandler &LocGlobMap() const override {
    LF_VERIFY_MSG(Mesh() != nullptr, "No valid FE space object: no mesh");
    LF_VERIFY_MSG(dofh_p_ != nullptr,
                  "No valid FE space object: no dof handler");
    return *dofh_p_;
  }

  /** @brief access to shape function layout for cells
   * @copydoc SclarFESpace::ShapeFunctionLayout(const lf::mesh::Entity&)
   */
  [[nodiscard]] std::shared_ptr<
      const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
  ShapeFunctionLayout(const lf::mesh::Entity &entity) const override;

  /** @brief access to shape function layout for cells
   *
   * @param ref_el_type type of entit, can be anything except for
   * lf::base::RefEl::kPoint()
   *
   * @warning NULL pointers may be returned by this method in case a finite
   * element specification was not given for a particular topological type of
   * entity.
   */
  [[nodiscard]] std::shared_ptr<
      const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
  ShapeFunctionLayout(lf::base::RefEl ref_el_type) const;

  /** @brief number of _interior_ shape functions associated to entities of
   * various types
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      const lf::mesh::Entity &entity) const override;

  /** @brief number of _interior_ shape functions associated to entities of
   * various types
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      lf::base::RefEl ref_el_type) const;

  /** @brief No special destructor */
  ~UniformScalarFESpace() override = default;

 private:
  /** Description of reference shape functions on triangular cells */
  std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
      rfs_tria_p_;
  /** Description of reference shape functions on quadrilateral cells */
  std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
      rfs_quad_p_;
  /** Description of reference shape functions on an edge */
  std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
      rfs_edge_p_;
  /** Description of refererence shape functions on a point */
  std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
      rfs_point_p_;
  /** Numbers of local shape functions for different types of entities */
  size_type num_rsf_node_{0}, num_rsf_edge_{0}, num_rsf_tria_{0},
      num_rsf_quad_{0};
  /** Local-to-global index map for the finite element space */
  std::unique_ptr<lf::assemble::UniformFEDofHandler> dofh_p_;

  /** Initialization of class member variables and consistency checks */
  void init(std::shared_ptr<const lf::mesh::Mesh> mesh_p);
  /** Checks whether some pointer are not valid */
  [[nodiscard]] bool check_ptr() const {
    LF_VERIFY_MSG(Mesh() != nullptr, "No valid FE space object: no mesh");
    LF_VERIFY_MSG(dofh_p_ != nullptr,
                  "No valid FE space object: no dof handler");
    LF_VERIFY_MSG((rfs_quad_p_ != nullptr) && (rfs_quad_p_ != nullptr),
                  "No valid FE space object: no rsfs for cells");
    return true;
  }

};  // end class definition UniformScalarFESpace

/** Output control variables for PrintInfo: */

/**
 * @brief mesh information will be printed
 *
 * To be used with
 * PrintInfo(std::ostream&, const UniformScalarFESpace<SCALAR>&, unsigned int)
 */
const unsigned int kUniformScalarFESpaceOutMesh = 1;

/**
 * @brief information about the dof handler will be printed.
 *
 * To be used with
 * PrintInfo(std::ostream&, const UniformScalarFESpace<SCALAR>&, unsigned int)
 */
const unsigned int kUniformScalarFESpaceOutDofh = 2;

/**
 * @brief information about the reference shape functions will be printed
 *
 * To be used with
 * PrintInfo(std::ostream&, const UniformScalarFESpace<SCALAR>&, unsigned int)
 */
const unsigned int kUniformScalarFESpaceOutRsfs = 4;

/**
 * @brief Print information about a UniformScalarFESpace to the given
 * stream object.
 * @param o The stream to which we should print
 * @param fes The UniformScalarFESpace that should be printed.
 * @param ctrl Controls the level of output:
 *
 * #### Level of output:
 * - if `ctrl & kUniformScalarFESpaceOutMesh`, mesh information will be printed.
 * - if `ctrl & kUniformScalarFESpaceOutDofh`, information about the dof handler
 * will be printed.
 * - if `ctrl & kUniformScalarFESpaceOutRsfs`, information about the reference
 * shape functions will be printed.
 *
 * @relates UniformScalarFESpace
 */
template <class SCALAR>
void PrintInfo(std::ostream &o, const UniformScalarFESpace<SCALAR> &fes,
               unsigned int ctrl = 0);

/**
 * @brief output operator for scalar parametric finite element space
 *
 * @relates UniformScalarFESpace
 */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o,
                         const UniformScalarFESpace<SCALAR> &fes);

// Initialization methods
template <typename SCALAR>
void UniformScalarFESpace<SCALAR>::init(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // Check validity and consistency of mesh pointer
  LF_VERIFY_MSG(mesh_p != nullptr, "Missing mesh!");
  LF_VERIFY_MSG((rfs_quad_p_ != nullptr) || (rfs_tria_p_ != nullptr),
                "Missing FE specification for cells");
  LF_VERIFY_MSG((mesh_p->DimMesh() == 2), "Only for 2D meshes");

  // Check whether all required finite element specifications are provided
  LF_VERIFY_MSG(
      (mesh_p->NumEntities(lf::base::RefEl::kTria()) == 0) ||
          (rfs_tria_p_ != nullptr),
      "Missing FE specification for triangles though mesh contains some");
  LF_VERIFY_MSG((mesh_p->NumEntities(lf::base::RefEl::kQuad()) == 0) ||
                    (rfs_quad_p_ != nullptr),
                "Missing FE specification for quads though mesh contains some");

  // Compatibility checks and initialization of numbers of shape functions
  // In particular only a single shape function may be associated to a node
  // in the case of a SCALAR finite element space
  if (rfs_tria_p_ != nullptr) {
    // Probe local shape functions on a triangle
    LF_VERIFY_MSG((*rfs_tria_p_).RefEl() == lf::base::RefEl::kTria(),
                  "Wrong type for triangle!");
    LF_VERIFY_MSG((*rfs_tria_p_).NumRefShapeFunctions(2) <= 1,
                  "At most one shape function can be assigned to each vertex");
    // Initialize numbers of shape functions associated to entities
    num_rsf_node_ = (*rfs_tria_p_).NumRefShapeFunctions(2);
    num_rsf_edge_ = (*rfs_tria_p_).NumRefShapeFunctions(1);
    num_rsf_tria_ = (*rfs_tria_p_).NumRefShapeFunctions(0);
  }
  if (rfs_quad_p_ != nullptr) {
    // Probe local shape functions for QUADs
    LF_VERIFY_MSG((*rfs_quad_p_).RefEl() == lf::base::RefEl::kQuad(),
                  "Wrong type for quad!");
    LF_VERIFY_MSG((*rfs_quad_p_).NumRefShapeFunctions(2) <= 1,
                  "At most one shape function can be assigned to each vertex");
    // Initialize numbers of shape functions associated to entities
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
  // edges. Those must be the same for all reference shape function descriptions
  // passed to the finite element space.
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
  dofh_p_ = std::make_unique<lf::assemble::UniformFEDofHandler>(
      std::move(mesh_p), rsf_layout);
}

template <typename SCALAR>
std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
UniformScalarFESpace<SCALAR>::ShapeFunctionLayout(
    const lf::mesh::Entity &entity) const {
  return ShapeFunctionLayout(entity.RefEl());
}

template <typename SCALAR>
std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
UniformScalarFESpace<SCALAR>::ShapeFunctionLayout(
    lf::base::RefEl ref_el_type) const {
  // Retrieve specification of local shape functions
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return rfs_point_p_;
    }
    case lf::base::RefEl::kSegment(): {
      // Null pointer valid return value: indicates that a shape function
      // description for edges is missing.
      // LF_ASSERT_MSG(rfs_edge_p_ != nullptr, "No RSF for edges!");
      return rfs_edge_p_;
    }
    case lf::base::RefEl::kTria(): {
      // Null pointer valid return value: indicates that a shape function
      // description for triangular cells is missing.
      // LF_ASSERT_MSG(rfs_tria_p_ != nullptr, "No RSF for triangles!");
      return rfs_tria_p_;
    }
    case lf::base::RefEl::kQuad(): {
      // Null pointer valid return value: indicates that a shape function
      // description for quadrilaterals is missing.
      // LF_ASSERT_MSG(rfs_quad_p_ != nullptr, "No RSF for quads!");
      return rfs_quad_p_;
    }
    default: {
      LF_VERIFY_MSG(false, "Illegal entity type");
    }
  }
  return nullptr;
}

template <typename SCALAR>
size_type UniformScalarFESpace<SCALAR>::NumRefShapeFunctions(
    const lf::mesh::Entity &entity) const {
  return NumRefShapeFunctions(entity.RefEl());
}

/* number of _interior_ shape functions associated to entities of various types
 */
template <typename SCALAR>
size_type UniformScalarFESpace<SCALAR>::NumRefShapeFunctions(
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

template <typename SCALAR>
void PrintInfo(std::ostream &o, const UniformScalarFESpace<SCALAR> &fes,
               unsigned ctrl) {
  o << "Uniform scalar FE space, dim = " << fes.LocGlobMap().NumDofs()
    << std::endl;

  if ((ctrl & kUniformScalarFESpaceOutMesh) != 0) {
    o << *fes.Mesh() << std::endl;
  }
  if ((ctrl & kUniformScalarFESpaceOutDofh) != 0) {
    o << fes.LocGlobMap() << std::endl;
  }
  if ((ctrl & kUniformScalarFESpaceOutRsfs) != 0) {
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kPoint()) << " rsfs @ nodes"
      << std::endl;
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kSegment())
      << " rsfs @ edges" << std::endl;
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kTria())
      << " rsfs @ triangles" << std::endl;
    o << fes.NumRefShapeFunctions(lf::base::RefEl::kQuad()) << " rsfs @ quads"
      << std::endl;
  }
}

/** output operator for scalar parametric finite element space */
template <typename SCALAR>
std::ostream &operator<<(std::ostream &o,
                         const UniformScalarFESpace<SCALAR> &fes) {
  fes.PrintInfo(o, 0);
  return o;
}

}  // namespace lf::uscalfe

#endif
