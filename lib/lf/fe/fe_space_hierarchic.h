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
 * @brief Lagrangian Finite Element Space of arbitrary degree
 *
 * @tparam SCALAR underlying scalar type, usually either `double` or
 * `complex<double>`
 *
 * This FE Space contains hierarchic Finite Elements, meaning that
 * the function spaces of lower polynomial degree are contained
 * in the higher order function spaces.
 *
 * The polynomial degree is constant over the whole mesh meaning
 * that there is no \f$p\f$-refinement supported.
 *
 * @note Some of the pointers may be NULL. For instance, if all computations
 *       are done on purely triangular meshes then a finite element
 * specification for quadrilaterals need not be given.
 */
template <typename SCALAR>
class FeSpaceHierarchic : public ScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  FeSpaceHierarchic() = delete;
  FeSpaceHierarchic(const FeSpaceHierarchic &) = delete;
  FeSpaceHierarchic(FeSpaceHierarchic &&) noexcept = default;
  FeSpaceHierarchic &operator=(const FeSpaceHierarchic &) = delete;
  FeSpaceHierarchic &operator=(FeSpaceHierarchic &&) noexcept = default;

  /**
   * @brief Constructor: Sets up the dof handler
   * @param mesh_p A shared pointer to the underlying mesh (immutable)
   * @param N The polynomial degree of the Finite Element Space
   */
  explicit FeSpaceHierarchic(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, unsigned N)
      : ScalarFESpace<SCALAR>(),
        mesh_p_(mesh_p),
        ref_el_(mesh_p),
        degree_(N),
        dofh_(Init()) {}

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
   * @copydoc ScalarFESpace::ShapeFunctionLayout(const lf::mesh::Entity&)
   */
  [[nodiscard]] ScalarReferenceFiniteElement<SCALAR> const *ShapeFunctionLayout(
      const lf::mesh::Entity &entity) const override {
    return std::visit(
        [](const auto &fe) -> ScalarReferenceFiniteElement<SCALAR> const * {
          if constexpr (std::is_same_v<decltype(fe), const std::monostate &>) {
            LF_VERIFY_MSG(false,
                          "Something is wrong, no "
                          "ScalarReferenceFiniteElement found for "
                          "this entity.");
            return nullptr;  // make compiler happy.
          } else {
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

  ~FeSpaceHierarchic() override = default;

 private:
  /* Underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;

  lf::mesh::utils::AllCodimMeshDataSet<
      std::variant<std::monostate, FePoint<SCALAR>, FeHierarchicSegment<SCALAR>,
                   FeHierarchicTria<SCALAR>, FeHierarchicQuad<SCALAR>>>
      ref_el_;
  unsigned degree_;
  lf::assemble::UniformFEDofHandler dofh_;

  // R.H. Detailed comment needed for this function
  [[nodiscard]] assemble::UniformFEDofHandler Init() {
    // Initialize all shape function layouts for nodes
    size_type num_rsf_node = 1;
    for (auto entity : Mesh()->Entities(2)) {
      ref_el_(*entity) = FePoint<SCALAR>(degree_);
    }
    // Initialize all shape function layouts for the edges
    size_type num_rsf_edge = 0;
    for (auto entity : Mesh()->Entities(1)) {
      FeHierarchicSegment<SCALAR> fe(degree_, entity->RelativeOrientations());
      num_rsf_edge = fe.NumRefShapeFunctions(0);
      ref_el_(*entity) = std::move(fe);
    }
    // Initialize all shape function layouts for the cells
    size_type num_rsf_tria = 0;
    size_type num_rsf_quad = 0;
    for (auto entity : Mesh()->Entities(0)) {
      switch (entity->RefEl()) {
        case lf::base::RefEl::kTria(): {
          FeHierarchicTria<SCALAR> fe(degree_, entity->RelativeOrientations());
          num_rsf_tria = fe.NumRefShapeFunctions(0);
          ref_el_(*entity) = std::move(fe);
          break;
        }
        case lf::base::RefEl::kQuad(): {
          FeHierarchicQuad<SCALAR> fe(degree_, entity->RelativeOrientations());
          num_rsf_quad = fe.NumRefShapeFunctions(0);
          ref_el_(*entity) = std::move(fe);
          break;
        }
        default:
          LF_VERIFY_MSG(false, "Illegal entity type");
      }
    }
    // Initialize the dofhandler
    lf::assemble::UniformFEDofHandler::dof_map_t rsf_layout{
        {lf::base::RefEl::kPoint(), num_rsf_node},
        {lf::base::RefEl::kSegment(), num_rsf_edge},
        {lf::base::RefEl::kTria(), num_rsf_tria},
        {lf::base::RefEl::kQuad(), num_rsf_quad}};
    return lf::assemble::UniformFEDofHandler(Mesh(), rsf_layout);
  }
};

}  // end namespace lf::fe

#endif  // LF_USCALFE_FE_SPACE_HP_H_
