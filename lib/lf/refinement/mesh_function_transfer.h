#ifndef LF_REFINEMENT_MESH_FUNCTION_TRANSFER_H
#define LF_REFINEMENT_MESH_FUNCTION_TRANSFER_H

#include <lf/assemble/dofhandler.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/uscalfe.h>
#include <Eigen/Dense>
#include <type_traits>

namespace lf::refinement {

/**
 * @brief A MeshFunction representing interpolation on a lf::refinement::MeshHierarchy
 * @tparam The type of mesh function to interpolate on a finer mesh
 */
template <typename MF>
class MeshFunctionTransfer {
  using mf_t = std::remove_cv_t<std::remove_reference_t<MF>>;
  static_assert(lf::mesh::utils::isMeshFunction<mf_t>,
                "MF is not a valid MeshFunction");

  // Metafunction to determine whether MF provides a `getMesh` method
  template <typename MF_TEST,
            typename = decltype(std::declval<MF_TEST>().getMesh())>
  static std::true_type has_getMesh_impl(int);
  template <typename MF_TEST>
  static std::false_type has_getMesh_impl(...);
  static constexpr bool provides_getMesh =
      decltype(has_getMesh_impl<mf_t>(int{})){};

 public:
  /**
   * @brief Copy Constructor
   */
  MeshFunctionTransfer(const MeshFunctionTransfer &) = default;

  /**
   * @brief Move Constructor
   */
  MeshFunctionTransfer(MeshFunctionTransfer &&) = default;

  /**
   * @brief Constructor
   * @param mh The lf::refinement::MeshHierarchy along which the function should be interpolated
   * @param mf The mesh function to be interpolated on a finer mesh
   * @param level_coarse The level of the Mesh which `mf` is defined on
   * @param level_fine The level of the Mesh this wrapper MeshFunction is defined on
   */
  MeshFunctionTransfer(const lf::refinement::MeshHierarchy &mh, const MF &mf,
                       lf::base::size_type level_coarse, lf::base::size_type level_fine)
      : mh_(mh), mf_(mf), level_coarse_(level_coarse), level_fine_(level_fine) {
    // Assert that the `level_coarse` parameter does not point to the finest
    // mesh
    LF_ASSERT_MSG(
        level_coarse < mh.NumLevels() - 1,
        "level_coarse must not point to the finest mesh in the hierarchy");
    // Assert that the `level_fine` parameter points to a finer mesh than `level_coarse`
    LF_ASSERT_MSG(level_fine > level_coarse, "level_fine must be bigger than level_fine");
    LF_ASSERT_MSG(level_fine < mh.NumLevels(), "level_fine must point to a valid mesh in the mesh hierarchy");
    // If the mesh function is defined over an FE space, assert the correctness
    // of the `level_coarse` parameter
    if constexpr (provides_getMesh) {
      LF_ASSERT_MSG(mh.getMesh(level_coarse) == mf.getMesh(),
                    "Invalid level_coarse provided");
    }
  }

  /**
   * @brief Evaluate the MeshFunction
   * @param e The Entity on which to evaluate this MeshFunction
   * @param local The local coordinates inside `e` where the MeshFunction should be evaluated
   * @returns A std::vector of point evaluations of the MeshFunction
   */
  decltype(auto) operator()(const lf::mesh::Entity &e,
                            const Eigen::MatrixXd &local) const {
    auto rel_geom = mh_.GeometryInParent(level_fine_, e);
    auto parent = mh_.ParentEntity(level_fine_, e);
    auto local_parent = rel_geom->Global(local);
    for (lf::base::size_type lvl = level_fine_-1 ; lvl > level_coarse_ ; --lvl) {
	rel_geom = mh_.GeometryInParent(lvl, *parent);
	parent = mh_.ParentEntity(lvl, *parent);
	local_parent = rel_geom->Global(local_parent);
    }
    return mf_(*parent, local_parent);
  }

  /**
   * @brief Access the underlying Mesh
   * @returns A shared pointer to the Mesh on which this MeshFunction is defined
   */
  std::shared_ptr<const lf::mesh::Mesh> getMesh() const {
    return mh_.getMesh(level_coarse_ + 1);
  }

 private:
  const lf::refinement::MeshHierarchy &mh_;
  const MF &mf_;
  const lf::base::size_type level_coarse_;
  const lf::base::size_type level_fine_;
};

template <typename MF>
MeshFunctionTransfer(const lf::refinement::MeshHierarchy &, const MF &,
                     lf::base::size_type, lf::base::size_type)
    ->MeshFunctionTransfer<MF>;

}  // namespace lf::refinement

#endif  // LF_REFINEMENT_MESH_FUNCTION_TRANSFER_H
