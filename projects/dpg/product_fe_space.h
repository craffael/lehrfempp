#ifndef PROJECTS_DPG_PRODUCT_FE_SPACE_H
#define PROJECTS_DPG_PRODUCT_FE_SPACE_H

/**
 * @file
 * @brief Data structure describing (cartesian)-product finite element spaces
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <vector>

#include <lf/base/base.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"
#include "product_dofhandler.h"

namespace projects::dpg {

/**
 * @headerfile projects/dpg/product_fe_space.h
 * @brief cartesian/prodcut space of finite element functions on a hybrid 2D
 * mesh.
 *
 * @tparam SCALAR underlying scalar type for all components, usually 'double'
 *
 *
 * This class provides functionality similar to the class
 *lf::uscalfe::UniformScalarFESpace for product/cartesian finite element spaces.
 *These spaces are of the form
 *
 *\f[ U = U_0 \times U_1 \times
 * \dots \times U_{n-1} \f] where each of the \f$ U_i \f$ is a  function space
 * of scalar  valued functions (see the description in
 *ProductUniformFEDofHandler).
 *
 * To allow this, this class provides some further member functions that give
 * extra information about the contained FE spaces. It furthermore uses a
 *ProductUniformFEDofHandler internally resulting in a different dof ordering.
 *
 * Furthermore, this class also weakens some constraints on the shape function
 *layouts compared to the lf::uscalfe::UniformScalarFESpace in order to
 *represent function spaces as they  occur in the description of DPG methods:
 *
 * - There may be multiple global shape functions associated with the same
 *vertex (usually they will belong to different components)
 * - Some components may not specify any local shape functions on triangles or
 * quadrilaterals, but only on segments.
 *
 *
 * @note Some of the pointers may be NULL, which indicates that the shape
 * function descriptions are missing.
 */
template <typename SCALAR>
class ProductUniformFESpace {
 public:
  using scalar = SCALAR;

  /** default constructor
   * @note creates an invalid object that can not be used.
   */
  ProductUniformFESpace() = default;
  ProductUniformFESpace(const ProductUniformFESpace &) = delete;
  ProductUniformFESpace(ProductUniformFESpace &&) noexcept = default;
  ProductUniformFESpace &operator=(const ProductUniformFESpace &) = delete;
  ProductUniformFESpace &operator=(ProductUniformFESpace &&) noexcept = default;

  /**
   * @brief Main construcotr: sets up the local-to-global index mapping
   * (dofhandler)
   * @param mesh_p shared pointer to underlying mesh (imutable)
   * @param rfs_tria_v vector of shared pointers to layout description of
   * reference shape functions on triangular cells for each component.
   * @param rfs_quad_v vector of shared pointers to layout description of
   * reference shape functions on quadrilateral cells for each component.
   * @param rfs_edge_v vector of shared pointers to layout description of
   * reference shape functions on the edges for each component.
   *
   * The schemes for local shape functions belonging to some component
   * have to satisfy the following compatibility conditions.
   *
   *  - The number of interior shape functions for edges of triangles and
   * quadrilaterals must agree.
   *
   * @note none of the shape function layouts needs to be specified; just pass
   * a null pointer. This will restrict the applicability
   * of the resulting finite element space object.
   */
  ProductUniformFESpace(
      std::shared_ptr<const lf::mesh::Mesh> mesh_p,
      std::vector<std::shared_ptr<
          const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
          rfs_tria_v,
      std::vector<std::shared_ptr<
          const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
          rfs_quad_v,
      std::vector<std::shared_ptr<
          const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
          rfs_edge_v,
      std::vector<std::shared_ptr<
          const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
          rfs_point_v)
      : mesh_p_(std::move(mesh_p)),
        rfs_tria_v_(std::move(rfs_tria_v)),
        rfs_quad_v_(std::move(rfs_quad_v)),
        rfs_edge_v_(std::move(rfs_edge_v)),
        rfs_point_v_(std::move(rfs_point_v)) {
    numComponents_ = rfs_tria_v_.size();
    init();
  }

  /** @brief acess to underlying mesh
   *  @return a shared _pointer_ to the mesh
   */
  [[nodiscard]] std::shared_ptr<const lf::mesh::Mesh> Mesh() const {
    LF_ASSERT_MSG(mesh_p_ != nullptr, "Invalid FE space, no mesh");
    return mesh_p_;
  }

  /** @brief access to associated local-to-global map
   * @return a reference to the ProductUniformFEDofHandler object (immutable)
   */
  [[nodiscard]] const ProductUniformFEDofHandler &LocGlobMap() const {
    LF_ASSERT_MSG(mesh_p_ != nullptr, "Invalid FE space, no mesh");
    LF_ASSERT_MSG(dofh_p_ != nullptr, "Invalid FE space, no dofhandler");
    return *dofh_p_;
  }

  /** @brief access to shape function layout
   *
   * @param ref_el_type type of entit, can be anything except for
   * lf::base::RefEl::kPoint()
   * @param component index of the component whose shape function
   * layout is returned.
   *
   * @warning NULL pointers may be returned by this method in case a finite
   * element specification was not given for a particular topological type of
   * entity and e particular component.
   */
  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>
  ShapeFunctionLayout(lf::base::RefEl ref_el_type, size_type component) const;

  /** @brief number of interior shape functions associated to entities of
   * various types for a given component */
  [[nodiscard]] size_type NumRefShapeFunctions(lf::base::RefEl ref_el_type,
                                               size_type component) const;

  /** @brief number of components of the product space */
  [[nodiscard]] size_type NumComponents() const { return numComponents_; }

  /** @brief Constructs a lf::uscalfe::UniformScalarFESpace that describes the
   * finite element space associated to a certain component
   * @note The resulting FESpace has no more information about the other
   * components
   * @warning For some components the construction of the
   * lf::uscalfe::UniformScalarFESpace may be impossible, since some of the
   * constraints regarding the shape function layouts were weakened.
   */
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<SCALAR>> ComponentFESpace(
      size_type component) {
    return std::make_shared<lf::uscalfe::UniformScalarFESpace<SCALAR>>(
        mesh_p_, rfs_tria_v_[component], rfs_quad_v_[component],
        rfs_edge_v_[component], rfs_point_v_[component]);
  }

  /** @brief no special destructor*/
  virtual ~ProductUniformFESpace() = default;

 private:
  /** underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;

  /** descritpions of reference shape functions for all components on
   *  different types of entities. */
  std::vector<
      std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_tria_v_;
  std::vector<
      std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_quad_v_;
  std::vector<
      std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_edge_v_;
  std::vector<
      std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_point_v_;

  /** numbers of local shape functions for all components
   * on different types of entities. */
  std::vector<size_type> num_rsf_node_;
  std::vector<size_type> num_rsf_edge_;
  std::vector<size_type> num_rsf_tria_;
  std::vector<size_type> num_rsf_quad_;

  /** local-to-global index map for the finite element space */
  std::unique_ptr<ProductUniformFEDofHandler> dofh_p_;

  /** number of components */
  size_type numComponents_ = 0;

  /** Initialization of class member variables and consistency checks */
  void init();
};

// Initialization method.
template <typename SCALAR>
void ProductUniformFESpace<SCALAR>::init() {
  // check validity and consistency of mesh pointer:
  LF_ASSERT_MSG(mesh_p_ != nullptr, "No mesh");
  LF_ASSERT_MSG(mesh_p_->DimMesh() == 2, "Only 2D cells");

  // we do not check, wheter a specification of reference finite elements
  // for cells (triangles and quadrilaterals) is given, since this may not
  // be desired necessary for some components in  a dpg method (fluxes)

  // check validity and consistency of the number of components:
  LF_ASSERT_MSG(numComponents_ != 0, "0 components FE space.");
  LF_ASSERT_MSG(numComponents_ == rfs_tria_v_.size() &&
                    numComponents_ == rfs_quad_v_.size() &&
                    numComponents_ == rfs_edge_v_.size() &&
                    numComponents_ == rfs_point_v_.size(),
                "Missmatch in number of shape functions passed for different "
                "entity tpyes.");

  // resize vectors saving the number of loal shape functions for each
  // component:
  num_rsf_tria_.resize(numComponents_);
  num_rsf_quad_.resize(numComponents_);
  num_rsf_edge_.resize(numComponents_);
  num_rsf_node_.resize(numComponents_);

  // compatibility checks and initialization of number of shape functions:
  for (size_type component = 0; component < numComponents_; component++) {
    num_rsf_tria_[component] = 0;
    num_rsf_quad_[component] = 0;
    num_rsf_edge_[component] = 0;
    num_rsf_node_[component] = 0;

    if (rfs_tria_v_[component] != nullptr) {
      // probe reference finite element
      LF_ASSERT_MSG(rfs_tria_v_[component]->RefEl() == lf::base::RefEl::kTria(),
                    "wrong type for triangle in component" << component);
      // initialize number of shape functions associated to entitites.
      num_rsf_node_[component] =
          rfs_tria_v_[component]->NumRefShapeFunctions(2);
      num_rsf_edge_[component] =
          rfs_tria_v_[component]->NumRefShapeFunctions(1);
      num_rsf_tria_[component] =
          rfs_tria_v_[component]->NumRefShapeFunctions(0);
    }
    if (rfs_quad_v_[component] != nullptr) {
      // probe reference finite element for quads.
      LF_ASSERT_MSG(rfs_quad_v_[component]->RefEl() == lf::base::RefEl::kQuad(),
                    "Wrong type for quad in component " << component);
      // initialize number of shape functions associated to entitites.
      num_rsf_node_[component] =
          rfs_quad_v_[component]->NumRefShapeFunctions(2);
      num_rsf_edge_[component] =
          rfs_quad_v_[component]->NumRefShapeFunctions(1);
      num_rsf_quad_[component] =
          rfs_quad_v_[component]->NumRefShapeFunctions(0);
    }
    if (rfs_edge_v_[component] != nullptr) {
      LF_ASSERT_MSG(
          rfs_edge_v_[component]->RefEl() == lf::base::RefEl::kSegment(),
          "Wrong type for edge in component " << component);
      num_rsf_node_[component] =
          rfs_edge_v_[component]->NumRefShapeFunctions(1);
      num_rsf_edge_[component] =
          rfs_edge_v_[component]->NumRefShapeFunctions(0);
    }

    // Compatibility check for numbers of local shape functions associated with
    // edges. Those must be the same for all reference shape function
    // descriptions in a component passed to the finite element space.
    if ((rfs_tria_v_[component] != nullptr) &&
        (rfs_quad_v_[component] != nullptr)) {
      LF_ASSERT_MSG((rfs_tria_v_[component]->NumRefShapeFunctions(2)) ==
                        (rfs_quad_v_[component]->NumRefShapeFunctions(2)),
                    "#RSF missmatch on nodes in component"
                        << component << ": "
                        << rfs_tria_v_[component]->NumRefShapeFunctions(2)
                        << "<->"
                        << rfs_quad_v_[component]->NumRefShapeFunctions(2));
      LF_ASSERT_MSG((rfs_tria_v_[component]->NumRefShapeFunctions(1)) ==
                        (rfs_quad_v_[component]->NumRefShapeFunctions(1)),
                    "#RSF missmatch on edges in component"
                        << component << ": "
                        << rfs_tria_v_[component]->NumRefShapeFunctions(1)
                        << "<->"
                        << rfs_quad_v_[component]->NumRefShapeFunctions(1));
    }
    if ((rfs_tria_v_[component] != nullptr) &&
        (rfs_edge_v_[component] != nullptr)) {
      LF_ASSERT_MSG((rfs_tria_v_[component]->NumRefShapeFunctions(2)) ==
                        (rfs_edge_v_[component]->NumRefShapeFunctions(1)),
                    "#RSF missmatch on nodes in component"
                        << component << ": "
                        << rfs_tria_v_[component]->NumRefShapeFunctions(2)
                        << "<->"
                        << rfs_edge_v_[component]->NumRefShapeFunctions(1));
      LF_ASSERT_MSG((rfs_tria_v_[component]->NumRefShapeFunctions(1)) ==
                        (rfs_edge_v_[component]->NumRefShapeFunctions(0)),
                    "#RSF missmatch on edges in component"
                        << component << ": "
                        << rfs_tria_v_[component]->NumRefShapeFunctions(1)
                        << "<->"
                        << rfs_edge_v_[component]->NumRefShapeFunctions(0));
    }
    if ((rfs_quad_v_[component] != nullptr) &&
        (rfs_edge_v_[component] != nullptr)) {
      LF_ASSERT_MSG((rfs_quad_v_[component]->NumRefShapeFunctions(2)) ==
                        (rfs_edge_v_[component]->NumRefShapeFunctions(1)),
                    "#RSF missmatch on nodes in component"
                        << component << ": "
                        << rfs_quad_v_[component]->NumRefShapeFunctions(2)
                        << "<->"
                        << rfs_edge_v_[component]->NumRefShapeFunctions(1));
      LF_ASSERT_MSG((rfs_quad_v_[component]->NumRefShapeFunctions(1)) ==
                        (rfs_edge_v_[component]->NumRefShapeFunctions(0)),
                    "#RSF missmatch on edges in component"
                        << component << ": "
                        << rfs_quad_v_[component]->NumRefShapeFunctions(1)
                        << "<->"
                        << rfs_edge_v_[component]->NumRefShapeFunctions(0));
    }
  }

  // initialization of dof handler
  // collect the number of interior reference shape functions for all
  // components
  std::vector<dof_map_t> rsf_layouts(numComponents_);
  for (size_type component = 0; component < numComponents_; component++) {
    rsf_layouts[component] = {
        {lf::base::RefEl::kPoint(), num_rsf_node_[component]},
        {lf::base::RefEl::kSegment(), num_rsf_edge_[component]},
        {lf::base::RefEl::kTria(), num_rsf_tria_[component]},
        {lf::base::RefEl::kQuad(), num_rsf_quad_[component]}};
  }
  dofh_p_ = std::make_unique<ProductUniformFEDofHandler>(mesh_p_, rsf_layouts);
}

template <typename SCALAR>
std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>
ProductUniformFESpace<SCALAR>::ShapeFunctionLayout(lf::base::RefEl ref_el_type,
                                                   size_type component) const {
  // Retrive specification of local shape functions for
  // a certain component.
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return rfs_point_v_[component];
    }
    case lf::base::RefEl::kSegment(): {
      // Null pointer valid return value:
      // indicates that a shape function description for edges is missing
      return rfs_edge_v_[component];
    }
    case lf::base::RefEl::kTria(): {
      // Null pointer valid return value:
      // indicates that a shape function description for triangles is missing.
      return rfs_tria_v_[component];
    }
    case lf::base::RefEl::kQuad(): {
      // Null pointer valid return value:
      // indicates that a shape function description for quadrilaterals is
      // missing.
      return rfs_quad_v_[component];
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal entity type");
    }
  }
  return nullptr;
}

template <typename SCALAR>
size_type ProductUniformFESpace<SCALAR>::NumRefShapeFunctions(
    lf::base::RefEl ref_el_type, size_type component) const {
  // Retrieve number of interior shape functions from rsf layouts
  // for a certain component.
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return num_rsf_node_[component];
    }
    case lf::base::RefEl::kSegment(): {
      return num_rsf_edge_[component];
    }
    case lf::base::RefEl::kTria(): {
      return num_rsf_tria_[component];
    }
    case lf::base::RefEl::kQuad(): {
      return num_rsf_quad_[component];
    }
    dafault : { LF_ASSERT_MSG(false, "Illegal entity type"); }
  }
  return 0;
}

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_PRODUCT_FE_SPACE_H
