#ifndef THESIS_POST_PROCESSING_MESH_FUNCTION_INTERPOLATION_H
#define THESIS_POST_PROCESSING_MESH_FUNCTION_INTERPOLATION_H

/**
 * @file mesh_function_interpolation.h
 * @brief Bring mesh functions down in a mesh hierarchy
 */

#include <lf/base/lf_exception.h>
#include <lf/refinement/mesh_hierarchy.h>

#include <vector>

namespace projects::ipdg_stokes {

namespace post_processing {

/**
 * @brief A MeshFunction representing another MeshFunction on a coarser mesh
 * @tparam INNER_MF The MeshFunction on the coarser mesh
 */
template <typename INNER_MF>
class MeshFunctionInterpolation {
 public:
  /**
   * @brief Create a new MeshFunction on refinement level `to` from a mesh
   * function on refinement level `from`
   * @param inner The MeshFunction on level `from`
   * @param meshes The MeshHierarchy the levels refer to
   * @param from The refinement level of the coarser mesh
   * @param to The ferinement level of the finer mesh
   *
   * Currently only meshes with affine geometries are supported. This could
   * change if a function lf::geometry::Geometry::Local(Eigen::MatrixXd &global)
   * would be provided.
   */
  MeshFunctionInterpolation(const INNER_MF &inner,
                            const lf::refinement::MeshHierarchy &meshes,
                            lf::base::size_type from, lf::base::size_type to)
      : inner_(inner),
        parents_(meshes.getMesh(to)->NumEntities(0)),
        mesh_to_(meshes.getMesh(to)) {
    // Initialize a map from entities in the fine mesh to entities in the coarse
    // mesh
    for (const auto cell : mesh_to_->Entities(0)) {
      parents_[mesh_to_->Index(*cell)] = parent(cell, meshes, from, to);
    }
  }

  decltype(auto) operator()(const lf::mesh::Entity &entity,
                            const Eigen::MatrixXd &local) {
    // Get the element in the coarser mesh which is the parent of entity
    const lf::mesh::Entity *parent = parents_[mesh_to_->Index(entity)];
    // Map the evaluation points from local coordinates in entity to local
    // coordinates in parent HACK: This assumes that the zeroth vertex of the
    // entity has local coordinates (0,0)
    const auto entity_geom = entity.Geometry();
    const auto parent_geom = parent->Geometry();
    if (!parent_geom->isAffine()) {
      throw lf::base::LfException(
          "MeshFunctionInterpolation only works with affine elements");
    }
    const Eigen::MatrixXd global = entity_geom->Global(local);
    const Eigen::MatrixXd global_shifted =
        global -
        parent_geom->Global(Eigen::MatrixXd::Zero(local.rows(), local.cols()));
    const Eigen::MatrixXd parent_local =
        parent_geom->JacobianInverseGramian(Eigen::Vector2d::Zero()) *
        global_shifted;
    // Evaluate the inner mesh function at the local coordinate inside parent
    return inner_(*parent, parent_local);
  }

 private:
  const INNER_MF inner_;
  std::vector<const lf::mesh::Entity *> parents_;
  std::shared_ptr<const lf::mesh::Mesh> mesh_to_;

  /**
   * @brief Find the parent entity given a child entity
   * @param entity A pointer to the cild entity
   * @param meshes The MeshHierarchy in which to find the parent
   * @param from The refinement level of the parent
   * @param to The refinement level of the cild entity
   */
  static const lf::mesh::Entity *parent(
      const lf::mesh::Entity *entity,
      const lf::refinement::MeshHierarchy &meshes, lf::base::size_type from,
      lf::base::size_type to) {
    if (from == to) {
      return entity;
    } else {
      const auto mesh_to = meshes.getMesh(to);
      const lf::mesh::Entity *immediate_parent =
          meshes.ParentInfos(to, 0)[mesh_to->Index(*entity)].parent_ptr;
      return parent(immediate_parent, meshes, from, to - 1);
    }
  }
};

}  // end namespace post_processing

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_POST_PROCESSING_MESH_FUNCTION_INTERPOLATION_H
