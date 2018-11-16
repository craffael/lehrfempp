/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief implementation of methods for data structures
 *        describing finite element spaces
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "fe_space.h"

namespace lf::fe {
// Constructor
UniformScalarFiniteElementSpace::UniformScalarFiniteElementSpace(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    std::unique_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::unique_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p)
    : mesh_p_(mesh_p),
      rfs_tria_p_(std::move(rfs_tria_p)),
      rfs_quad_p_(std::move(rfs_quad_p)) {
  LF_ASSERT_MSG((mesh_p_->DimMesh() == 2) && (rfs_tria_p_->Dimension() == 2) &&
                    (rfs_quad_p_->Dimension() == 2),
                "Only for 2D meshes");
  // Compatibility check for numbers of local shape functions associated with
  // sub-entitities
  LF_ASSERT_MSG(((*rfs_tria_p_).RefEl() == lf::base::RefEl::kTria()) ||
                    ((*rfs_quad_p_).RefEl() == lf::base::RefEl::kQuad()),
                "Unexpected type of reference cell");
  LF_ASSERT_MSG(((*rfs_tria_p_).NumRefShapeFunctions(2) == 1) &&
                    ((*rfs_quad_p_).NumRefShapeFunctions(2) == 1),
                "Exactly one shape function must be assigned to each vertex");
  LF_ASSERT_MSG(((*rfs_tria_p_).NumRefShapeFunctions(1) ==
                 (*rfs_quad_p_).NumRefShapeFunctions(1)),
                "#RSF mismatch on edges "
                    << (*rfs_tria_p_).NumRefShapeFunctions(1) << " <-> "
                    << (*rfs_quad_p_).NumRefShapeFunctions(1));
  // Initialization of dof handler starting with collecting the number of
  // interior reference shape functions
  lf::assemble::UniformFEDofHandler::dof_map_t rsf_layout{
      {lf::base::RefEl::kPoint(),
       NumRefShapeFunctions(lf::base::RefEl::kPoint())},
      {lf::base::RefEl::kSegment(),
       NumRefShapeFunctions(lf::base::RefEl::kSegment())},
      {lf::base::RefEl::kTria(),
       NumRefShapeFunctions(lf::base::RefEl::kTria())},
      {lf::base::RefEl::kQuad(),
       NumRefShapeFunctions(lf::base::RefEl::kQuad())}};
  dofh_p_ =
      std::make_unique<lf::assemble::UniformFEDofHandler>(mesh_p_, rsf_layout);
}

std::pair<const ScalarReferenceFiniteElement<double> &,
          const ScalarReferenceFiniteElement<double> &>
UniformScalarFiniteElementSpace::ShapeFunctionLayouts() const {
  LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
  LF_VERIFY_MSG(dofh_p_ != nullptr, "No valid FE space object: no dof handler");
  LF_VERIFY_MSG((rfs_quad_p_ != nullptr) && (rfs_quad_p_ != nullptr),
                "No valid FE space object: no rsfs");
  return {*rfs_tria_p_, *rfs_quad_p_};
}
/* access to shape function layout for triangular cells */
const ScalarReferenceFiniteElement<double>
    &UniformScalarFiniteElementSpace::TriaShapeFunctionLayout() const {
  check_ptr();
  return *rfs_tria_p_;
}
/* access to shape function layout for quadrilateral cells */
const ScalarReferenceFiniteElement<double>
    &UniformScalarFiniteElementSpace::QuadShapeFunctionLayout() const {
  check_ptr();
  return *rfs_quad_p_;
}

/* number of _interior_ shape functions associated to entities of various types
 */
size_type UniformScalarFiniteElementSpace::NumRefShapeFunctions(
    lf::base::RefEl ref_el_type) const {
  LF_VERIFY_MSG((rfs_quad_p_ != nullptr) && (rfs_quad_p_ != nullptr),
                "No valid FE space object: no rsfs");
  // Retrieve number of interior shape functions from rsf layouts
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return 1;
    }
    case lf::base::RefEl::kSegment(): {
      return (rfs_tria_p_->NumRefShapeFunctions(1));
    }
    case lf::base::RefEl::kTria(): {
      return (rfs_tria_p_->NumRefShapeFunctions(0));
    }
    case lf::base::RefEl::kQuad(): {
      return (rfs_quad_p_->NumRefShapeFunctions(0));
    }
    dafault : { LF_VERIFY_MSG(false, "Illegal entity type"); }
  }
  return 0;
}

/** output operator for scalar parametric finite element space */
std::ostream &operator<<(std::ostream &o,
                         const UniformScalarFiniteElementSpace &fes) {
  o << "Uniform scalar FE space, dim = " << fes.LocGlobMap().NoDofs()
    << std::endl;
  if (UniformScalarFiniteElementSpace::ctrl_ &
      UniformScalarFiniteElementSpace::kout_mesh) {
    o << fes.Mesh() << std::endl;
  }
  if (UniformScalarFiniteElementSpace::ctrl_ &
      UniformScalarFiniteElementSpace::kout_dofh) {
    o << fes.LocGlobMap() << std::endl;
  }
  if (UniformScalarFiniteElementSpace::ctrl_ &
      UniformScalarFiniteElementSpace::kout_rsfs) {
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
