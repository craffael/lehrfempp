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
// Output control variable
CONTROLDECLARECOMMENT(
    UniformScalarFiniteElementSpace, ctrl_,
    "UniformScalarFiniteElementSpace_ctrl",
    "Output control for class UniformScalarFiniteElementSpace");

// Initialization methods
void UniformScalarFiniteElementSpace::init() {
  LF_VERIFY_MSG((mesh_p_->DimMesh() == 2), "Only for 2D meshes");
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

const ScalarReferenceFiniteElement<double>
    &UniformScalarFiniteElementSpace::ShapeFunctionLayout(
        lf::base::RefEl ref_el_type) const {
  // Retrieve specification of local shape functions
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      LF_VERIFY_MSG(false, "No shape functions for point!");
      break;
    }
    case lf::base::RefEl::kSegment(): {
      LF_ASSERT_MSG(rfs_edge_p_ != nullptr, "No RSF for edges!");
      return *rfs_edge_p_;
    }
    case lf::base::RefEl::kTria(): {
      LF_ASSERT_MSG(rfs_tria_p_ != nullptr, "No RSF for triangles!");
      return *rfs_tria_p_;
    }
    case lf::base::RefEl::kQuad(): {
      LF_ASSERT_MSG(rfs_quad_p_ != nullptr, "No RSF for quads!");
      return *rfs_quad_p_;
    }
    dafault : { LF_VERIFY_MSG(false, "Illegal entity type"); }
  }
  return (*rfs_tria_p_);
}

/* number of _interior_ shape functions associated to entities of various types
 */
size_type UniformScalarFiniteElementSpace::NumRefShapeFunctions(
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
