/**
 * @file
 * @brief Lagrangian finite elements of arbitrary polynomials degree
 * @author Tobias Rohner
 * @date   May 20202
 * @copyright MIT License
 */

#ifndef LF_USCALFE_FE_SPACE_HP_H_
#define LF_USCALFE_FE_SPACE_HP_H_

#include <lf/mesh/utils/all_codim_mesh_data_set.h>

#include "hierarchic_fe.h"
#include "scalar_fe_space.h"

namespace lf::fe {

/**
 * @headerfile lf/fe/fe.h
 * @brief Finite Element Space that supports arbitrary, local degrees.
 *
 * @tparam SCALAR underlying scalar type, usually either `double` or
 * `complex<double>`
 *
 * This FE Space contains hierarchic Finite Elements, meaning that
 * the function spaces of lower polynomial degree are contained
 * in the higher order function spaces.
 *
 * The polynomial degree can vary from entity to entity, i.e. local
 * \f$p\f$-refinement is supported.
 *
 * A complete description of the basis functions and their dual basis can be
 * found <a
 * href="https://raw.githubusercontent.com/craffael/lehrfempp/master/doc/pfem/hierarchical_basis.pdf"
 * target="_blank"><b>here</b></a>.
 *
 * ### Example usage
 * The following code snippet computes the solution of the BVP
 * \f{align}
 * - \Delta u &= 1 && \text{on }\Omega := [0,1]^2 \\
 * u &= 0 && \text{on }\partial \Omega
 * \f}
 *
 * @snippet hierarchic_scalar_fe_space_snippets.cc Laplace
 */
template <typename SCALAR>
class HierarchicScalarFESpace : public ScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  HierarchicScalarFESpace() = delete;
  HierarchicScalarFESpace(const HierarchicScalarFESpace &) = delete;
  HierarchicScalarFESpace(HierarchicScalarFESpace &&) noexcept = default;
  HierarchicScalarFESpace &operator=(const HierarchicScalarFESpace &) = delete;
  HierarchicScalarFESpace &operator=(HierarchicScalarFESpace &&) noexcept =
      default;

  /**
   * @brief Construct a new Hierarchic FESpace with uniform polynomial degree.
   * @param mesh_p A shared pointer to the underlying mesh (immutable)
   * @param degree The uniform polynomial degree.
   */
  explicit HierarchicScalarFESpace(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, unsigned degree)
      : HierarchicScalarFESpace(mesh_p,
                                [degree](const mesh::Entity &e) -> unsigned {
                                  if (e.RefEl() == base::RefEl::kPoint()) {
                                    return 1;
                                  }
                                  return degree;
                                }) {}

  /**
   * @brief Construct a new Hierarchic FESpace with (possibly) varying
   * polynomial degrees.
   * @tparam F type of the `degree_functor`
   * @param mesh_p A shared pointer to the underlying mesh (immutable)
   * @param degree_functor A function object that assigns a polynomial degree to
   * every entity in the mesh. See below for more info.
   *
   * ### Degree Functor
   * The `degree_functor` must overload the call operator as follows:
   * ```
   * unsigned degree_functor(const lf::mesh::Entity& e)
   * ```
   * and should return the polynomial degree for the respective entity.
   * The degree functor will be called for all entities (i.e. edges, triangles,
   * quadrilaterals) in the mesh except for points.
   */
  template <class F, class = std::enable_if_t<
                         std::is_invocable_v<F, const mesh::Entity &>>>
  explicit HierarchicScalarFESpace(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, F &&degree_functor)
      : ScalarFESpace<SCALAR>(),
        mesh_p_(mesh_p),
        ref_el_(mesh_p),
        dofh_(mesh_p,
              [&degree_functor](const mesh::Entity &e) -> base::size_type {
                if (e.RefEl() == base::RefEl::kPoint()) {
                  return 1;
                }
                auto degree = degree_functor(e);
                switch (e.RefEl()) {
                  case base::RefEl::kSegment():
                    return degree - 1;
                  case base::RefEl::kTria():
                    return degree <= 2 ? 0 : (degree - 2) * (degree - 1) / 2;
                  case base::RefEl::kQuad():
                    return (degree - 1) * (degree - 1);
                  default:
                    LF_VERIFY_MSG(false, "Something went wrong.");
                }
              }) {
    Init(std::forward<F>(degree_functor));
  }

  /** @brief access to underlying mesh
   *  @return a shared _pointer_ to the mesh
   */
  [[nodiscard]] std::shared_ptr<const lf::mesh::Mesh> Mesh() const override {
    LF_VERIFY_MSG(mesh_p_ != nullptr, "No valid FE space object: no mesh");
    return mesh_p_;
  }

  /** @brief access to associated local-to-global map
   * @return a reference to the lf::assemble::DofHandler object (immutable)
   */
  [[nodiscard]] const lf::assemble::DofHandler &LocGlobMap() const override {
    LF_VERIFY_MSG(Mesh() != nullptr, "No valid FE space object: no mesh");
    return dofh_;
  }

  /** @brief access to shape function layout for cells
   * @copydoc ScalarFESpace::ShapeFunctionLayout()
   */
  [[nodiscard]] ScalarReferenceFiniteElement<SCALAR> const *ShapeFunctionLayout(
      const lf::mesh::Entity &entity) const override {
    return std::visit(
        [](const auto &fe) -> ScalarReferenceFiniteElement<SCALAR> const * {
          if constexpr (std::is_same_v<decltype(fe),
                                       const std::monostate &>) {  // NOLINT
            LF_VERIFY_MSG(false,
                          "Something is wrong, no "
                          "ScalarReferenceFiniteElement found for "
                          "this entity.");
            return nullptr;  // make compiler happy.
          } else {           // NOLINT
            return &fe;
          }
        },
        ref_el_(entity));
  }

  /** @brief number of _interior_ shape functions associated to entities of
   * various types
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      const lf::mesh::Entity &entity) const override {
    return ShapeFunctionLayout(entity)->NumRefShapeFunctions();
  }

  ~HierarchicScalarFESpace() override = default;

 private:
  /* Underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  quad::QuadRuleCache qr_cache_;
  // Store the ScalarReferenceFiniteElement for every entity.
  lf::mesh::utils::AllCodimMeshDataSet<
      std::variant<std::monostate, FePoint<SCALAR>, FeHierarchicSegment<SCALAR>,
                   FeHierarchicTria<SCALAR>, FeHierarchicQuad<SCALAR>>>
      ref_el_;
  lf::assemble::DynamicFEDofHandler dofh_;

  // R.H. Detailed comment needed for this function
  template <class F>
  void Init(F &&degree_functor) {
    // Initialize all shape function layouts for nodes
    const size_type num_rsf_node = 1;
    for (const auto *entity : mesh_p_->Entities(2)) {
      ref_el_(*entity) = FePoint<SCALAR>(1);
    }
    // Initialize all shape function layouts for the edges
    for (const auto *entity : mesh_p_->Entities(1)) {
      FeHierarchicSegment<SCALAR> fe(degree_functor(*entity), qr_cache_);

      ref_el_(*entity) = std::move(fe);
    }
    // Initialize all shape function layouts for the cells
    for (const auto *entity : mesh_p_->Entities(0)) {
      switch (entity->RefEl()) {
        case lf::base::RefEl::kTria(): {
          const std::array<unsigned, 3> edge_degrees{
              {degree_functor(*entity->SubEntities(1)[0]),
               degree_functor(*entity->SubEntities(1)[1]),
               degree_functor(*entity->SubEntities(1)[2])}};
          FeHierarchicTria<SCALAR> fe(degree_functor(*entity), edge_degrees,
                                      qr_cache_,
                                      entity->RelativeOrientations());
          ref_el_(*entity) = std::move(fe);
          break;
        }
        case lf::base::RefEl::kQuad(): {
          const std::array<unsigned, 4> edge_degrees{
              {degree_functor(*entity->SubEntities(1)[0]),
               degree_functor(*entity->SubEntities(1)[1]),
               degree_functor(*entity->SubEntities(1)[2]),
               degree_functor(*entity->SubEntities(1)[3])}};
          FeHierarchicQuad<SCALAR> fe(degree_functor(*entity), edge_degrees,
                                      qr_cache_,
                                      entity->RelativeOrientations());
          ref_el_(*entity) = std::move(fe);
          break;
        }
        default:
          LF_VERIFY_MSG(false, "Illegal entity type");
      }
    }
  }
};

}  // end namespace lf::fe

#endif  // LF_USCALFE_FE_SPACE_HP_H_
