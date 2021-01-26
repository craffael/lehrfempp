#ifndef PROJECTS_DPG_PRODUCT_FE_SPACE_FACTORY
#define PROJECTS_DPG_PRODUCT_FE_SPACE_FACTORY

/**
 * @file
 * @brief Factory class to create ProductUniformFEspaces.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include "discontinuous_fe_constant.h"
#include "discontinuous_scalar_reference_finite_element.h"
#include "dpg.h"
#include "product_fe_space.h"
#include "trace_scalar_reference_finite_element.h"

namespace projects::dpg {

/**
 * @brief Factory class to build  a ProductUniformFESpace.
 *
 * This class can be used to construct a  ProductUniformFESpace. The class
 * provides several methods to add new components  based on different continuity
 * properties as well as the degree of the local shape functiosn .
 */
template <typename SCALAR>
class ProductUniformFESpaceFactory {
 public:
  /** @brief standard constructors */
  ProductUniformFESpaceFactory(const ProductUniformFESpaceFactory&) = delete;
  ProductUniformFESpaceFactory(ProductUniformFESpaceFactory&&) noexcept =
      delete;
  ProductUniformFESpaceFactory& operator=(const ProductUniformFESpaceFactory&) =
      delete;
  ProductUniformFESpaceFactory& operator=(ProductUniformFESpaceFactory&&) =
      delete;

  /**
   * @brief main constructor, constructs a new builder
   * @param mesh_p the underlying mesh, on which the ProductUniformFESpace will
   * be built.
   */
  explicit ProductUniformFESpaceFactory(
      std::shared_ptr<const lf::mesh::Mesh> mesh_p)
      : mesh_p_(std::move(mesh_p)) {}

  /**
   * @brief Adds a \f$H^1(\Omega) \f$-conforming finite element to the space.
   * @param degree the polynomial degree for the local shape functions
   * @return The zero based index of the added component.
   *
   * Since \f$ H^1(\Omega) \f$-conforming finite elements are well
   * defined on mesh interfaces, reference shape functions on triangles,
   * quadrilaterals and segments will be specified.
   *
   */
  size_type AddH1Component(size_type degree);

  /**
   * @brief Adds a \f$L^2(\Omega) \f$-conforming finite element to the space.
   * @param degree the polynomial degree for the local shape functions
   * @return The zero based index of the added component.
   *
   * Since \f$L^2(\Omega) \f$-conforming finite elements are no longer well
   * defined on mesh interfaces, reference shape functions will only be
   * specified on triangles and quadrilaterals and no longer on segments.
   *
   * @note Since the local shape functions are polynomials this function can
   * also be used for \f$ H^1(\mathcal M) \f$ or \f$ H(\text{div}, \mathcal M)
   * \f$-conforming finite elements.
   */
  size_type AddL2Component(size_type degree);

  /**
   * @brief Adds a \f$H^{-1/2}(S) \f$-conforming finite element to the space,
   * that  represents a flux on the mesh skeleton
   * @param degree the polynomial degree for the local shape functions
   * @return The zero based index of the added component.
   *

   * Since fluxes only live on the mesh skeleton, local
   * reference shape functions are only specified on  segments and not  on
   triangles or
   * quadrilaterals.
   */
  size_type AddFluxComponent(size_type degree);

  /**
   * @brief Adds a \f$H^{1/2}(S) \f$-conforming finite element to the space,
   * that  represents a trace on the mesh skeleton.
   * @param degree the polynomial degree for the local shape functions
   * @return The zero based index of the added component.
   *
   * Since
   * traces live on the mesh skeleton,  local reference shape
   * functions on segments are specified. To to simplify assembly procedures
   * local shape functions on  triangles and quadrilaterals using the
   * TraceScalarReferenceFiniteElement wrapper are also provided.
   */
  size_type AddTraceComponent(size_type degree);

  /**
   * @brief Build the ProductUniformFespace based on the provided specifications

   */
  std::shared_ptr<ProductUniformFESpace<SCALAR>> Build();

  /**
   * @brief return the number of components that have already been added.
   */
  size_type NumComponents() { return num_components_; }

  virtual ~ProductUniformFESpaceFactory() = default;

 private:
  /**
   * @brief returns Lagrangian shape functions of a specified degree on a
   * specific reference element.
   * @param ref_el the reference element for which the shape functions should be
   * defined
   * @param degree the degree of the shape functions
   * @return returns a Pointer to an object representing Lagrangian shape
   * functions of the specified degree and on the specified reference element.
   */
  std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
  ReferenceFiniteElement(lf::base::RefEl ref_el, size_type degree);

  /** the underlying mesh on which the fe-space will be built. */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;

  /** descriptions of reference shape functions for already added components on
   * different types of entities */
  std::vector<
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_tria_v_;
  std::vector<
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_quad_v_;
  std::vector<
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_segment_v_;
  std::vector<
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>>
      rfs_point_v_;

  /** the number of already added components. */
  size_type num_components_ = 0;
};

template <typename SCALAR>
std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>>
ProductUniformFESpaceFactory<SCALAR>::ReferenceFiniteElement(
    lf::base::RefEl ref_el, size_type degree) {
  switch (ref_el) {
    case (lf::base::RefEl::kTria()): {
      switch (degree) {
        case (0):
          return std::make_shared<FeDiscontinuousO0Tria<SCALAR>>();
        case (1):
          return std::make_shared<lf::uscalfe::FeLagrangeO1Tria<SCALAR>>();
        case (2):
          return std::make_shared<lf::uscalfe::FeLagrangeO2Tria<SCALAR>>();
        case (3):
          return std::make_shared<lf::uscalfe::FeLagrangeO3Tria<SCALAR>>();
        default:
          LF_ASSERT_MSG(false, "no fe of degree " << degree
                                                  << " available on ref_el "
                                                  << ref_el);
      }
    }
    case (lf::base::RefEl::kQuad()): {
      switch (degree) {
        case (0):
          return std::make_shared<FeDiscontinuousO0Quad<SCALAR>>();
        case (1):
          return std::make_shared<lf::uscalfe::FeLagrangeO1Quad<SCALAR>>();
        case (2):
          return std::make_shared<lf::uscalfe::FeLagrangeO2Quad<SCALAR>>();
        case (3):
          return std::make_shared<lf::uscalfe::FeLagrangeO3Quad<SCALAR>>();
        default:
          LF_ASSERT_MSG(false, "no fe of degree " << degree
                                                  << "available on ref_el "
                                                  << ref_el);
      }
    }
    case (lf::base::RefEl::kSegment()): {
      switch (degree) {
        case (0):
          return std::make_shared<FeDiscontinuousO0Segment<SCALAR>>();
        case (1):
          return std::make_shared<lf::uscalfe::FeLagrangeO1Segment<SCALAR>>();
        case (2):
          return std::make_shared<lf::uscalfe::FeLagrangeO2Segment<SCALAR>>();
        case (3):
          return std::make_shared<lf::uscalfe::FeLagrangeO3Segment<SCALAR>>();
        default:
          LF_ASSERT_MSG(false, "no fe of degree " << degree
                                                  << "available on ref_el "
                                                  << ref_el);
      }
    }
    case (lf::base::RefEl::kPoint()): {
      return std::make_shared<lf::fe::FePoint<SCALAR>>(degree);
    }
    default: {
      LF_ASSERT_MSG(false, "unsupported reference element " << ref_el);
    }
  }
}

template <typename SCALAR>
size_type ProductUniformFESpaceFactory<SCALAR>::AddH1Component(
    size_type degree) {
  LF_ASSERT_MSG(
      degree > 0,
      "Can only represent H1 conforming finite element with degree > 0");
  // add continuous reference shape functions on triangles, quads and segments.
  rfs_tria_v_.push_back(
      ReferenceFiniteElement(lf::base::RefEl::kTria(), degree));
  rfs_quad_v_.push_back(
      ReferenceFiniteElement(lf::base::RefEl::kQuad(), degree));
  rfs_segment_v_.push_back(
      ReferenceFiniteElement(lf::base::RefEl::kSegment(), degree));
  rfs_point_v_.push_back(nullptr);
  return num_components_++;
}

template <typename SCALAR>
size_type ProductUniformFESpaceFactory<SCALAR>::AddL2Component(
    size_type degree) {
  // add discontinuous shape functions to triangles and quads.
  rfs_tria_v_.push_back(
      std::make_shared<DiscontinuousScalarReferenceFiniteElement<SCALAR>>(
          ReferenceFiniteElement(lf::base::RefEl::kTria(), degree)));
  rfs_quad_v_.push_back(
      std::make_shared<DiscontinuousScalarReferenceFiniteElement<SCALAR>>(
          ReferenceFiniteElement(lf::base::RefEl::kQuad(), degree)));
  rfs_segment_v_.push_back(nullptr);
  rfs_point_v_.push_back(nullptr);
  return num_components_++;
}

template <typename SCALAR>
size_type ProductUniformFESpaceFactory<SCALAR>::AddFluxComponent(
    size_type degree) {
  // add discontinuous shape functions to segments
  rfs_tria_v_.push_back(nullptr);
  rfs_quad_v_.push_back(nullptr);
  rfs_segment_v_.push_back(
      std::make_shared<DiscontinuousScalarReferenceFiniteElement<SCALAR>>(
          ReferenceFiniteElement(lf::base::RefEl::kSegment(), degree)));
  rfs_point_v_.push_back(nullptr);
  return num_components_++;
}

template <typename SCALAR>
size_type ProductUniformFESpaceFactory<SCALAR>::AddTraceComponent(
    size_type degree) {
  // add trace representing shape functions to triangles and quadrilaterals.
  // add full shape function spec. to segments.
  rfs_tria_v_.push_back(
      std::make_shared<TraceScalarReferenceFiniteElement<SCALAR>>(
          ReferenceFiniteElement(lf::base::RefEl::kTria(), degree)));
  rfs_quad_v_.push_back(
      std::make_shared<TraceScalarReferenceFiniteElement<SCALAR>>(
          ReferenceFiniteElement(lf::base::RefEl::kQuad(), degree)));
  rfs_segment_v_.push_back(
      ReferenceFiniteElement(lf::base::RefEl::kSegment(), degree));
  rfs_point_v_.push_back(nullptr);
  return num_components_++;
}

template <typename SCALAR>
std::shared_ptr<ProductUniformFESpace<SCALAR>>
ProductUniformFESpaceFactory<SCALAR>::Build() {
  // construct the fe space based on the provided local reference shape function
  // layouts.
  auto fe_space = std::make_shared<ProductUniformFESpace<SCALAR>>(
      mesh_p_, rfs_tria_v_, rfs_quad_v_, rfs_segment_v_, rfs_point_v_);

  // clear supplied information:
  rfs_tria_v_ = {};
  rfs_quad_v_ = {};
  rfs_segment_v_ = {};
  rfs_point_v_ = {};
  return fe_space;
}

}  // namespace projects::dpg
#endif  // PROJECTS_DPG_PRODUCT_FE_SPACE_FACTORY
