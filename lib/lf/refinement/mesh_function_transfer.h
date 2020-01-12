#ifndef LF_REFINEMENT_MESH_FUNCTION_TRANSFER_H
#define LF_REFINEMENT_MESH_FUNCTION_TRANSFER_H

#include <lf/assemble/dofhandler.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/mesh/utils/utils.h>
#include <Eigen/Dense>
#include <type_traits>

namespace lf::refinement {

template <typename MF>
class MeshFunctionTransfer {
  using mf_t = std::remove_cv_t<std::remove_reference_t<MF>>;
  static_assert(lf::mesh::utils::isMeshFunction<mf_t>,
                "MF is not a valid MeshFunction");

  // Metafunction to determine whether MF provides a `getMesh` method
  template <typename MF_TEST,
            typename = decltype(std::declval<MF_TEST>().getMesh())>
  static std::true_type has_getMesh_impl();
  template <typename MF_TEST>
  static std::false_type has_getMesh_impl();
  static constexpr bool provides_getMesh = decltype(has_getMesh_impl<mf_t>()){};

 public:
  MeshFunctionTransfer(const MeshFunctionTransfer &) = default;
  MeshFunctionTransfer(MeshFunctionTransfer &&) = default;
  MeshFunctionTransfer(const lf::refinement::MeshHierarchy &mh, const MF &mf,
                       lf::base::size_type level_coarse)
      : mh_(mh), mf_(mf), level_coarse_(level_coarse) {
    // Assert that the `level_coarse` parameter does not point to the finest
    // mesh
    LF_ASSERT_MSG(
        level_coarse < mh.NumLevels() - 1,
        "level_coarse must not point to the finest mesh in the hierarchy");
    // If the mesh function is defined over an FE space, assert the correctness
    // of the `level_coarse` parameter
    if constexpr (provides_getMesh) {
      LF_ASSERT_MSG(mh.getMesh(level_coarse) == mf.getMesh(),
                    "Invalid level_coarse provided");
    }
  }

  decltype(auto) operator()(const lf::mesh::Entity &e,
                            const Eigen::MatrixXd &local) const {
    const auto parent = mh_.ParentEntity(level_coarse_ + 1, e);
    const auto rel_geom = mh_.GeometryInParent(level_coarse_ + 1, e);
    const auto local_parent = rel_geom->Global(local);
    return mf_(parent, local_parent);
  }

  decltype(auto) getMesh() const { return mh_.getMesh(level_coarse_ + 1); }

 private:
  const lf::refinement::MeshHierarchy &mh_;
  const MF &mf_;
  const lf::base::size_type level_coarse_;
};

template <typename MF>
MeshFunctionTransfer(const lf::refinement::MeshHierarchy &, const MF &,
                     lf::base::size_type)
    ->MeshFunctionTransfer<MF>;

}  // namespace lf::refinement

#endif  // LF_REFINEMENT_MESH_FUNCTION_TRANSFER_H
