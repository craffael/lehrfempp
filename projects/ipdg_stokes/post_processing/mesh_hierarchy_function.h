#ifndef THESIS_POST_PROCESSING_MESH_HIERARCHY_FUNCTION_H
#define THESIS_POST_PROCESSING_MESH_HIERARCHY_FUNCTION_H

#include <functional>
#include <map>
#include <numeric>

#include <lf/refinement/mesh_hierarchy.h>

namespace projects::ipdg_stokes {

namespace post_processing {

template <typename Ret>
std::map<lf::base::size_type,
         std::function<Ret(const lf::mesh::Entity &, const Eigen::Vector2d &)>>
bringToFinestMesh(
    const lf::refinement::MeshHierarchy &meshes,
    const std::map<lf::base::size_type,
                   std::function<Ret(const lf::mesh::Entity &,
                                     const Eigen::Vector2d &)>> &f) {
  // Build a data structure storing the parent information for each entity in
  // the mesh hierarchy
  const lf::base::size_type num_cells =
      meshes.getMesh(meshes.NumLevels() - 1)->NumEntities(0);
  std::vector<std::vector<lf::base::size_type>> is_part_of(
      meshes.NumLevels(), std::vector<lf::base::size_type>(num_cells));
  std::iota(is_part_of[meshes.NumLevels() - 1].begin(),
            is_part_of[meshes.NumLevels() - 1].end(), 0);
  for (lf::base::size_type lvl = meshes.NumLevels() - 1; lvl > 0; --lvl) {
    const auto mesh = meshes.getMesh(lvl);
    const auto &child_info = meshes.CellChildInfos(lvl - 1);
    std::vector<lf::base::size_type> parent_info(mesh->NumEntities(0));
    for (lf::base::size_type i = 0; i < child_info.size(); ++i)
      for (auto child : child_info[i].child_cell_idx) parent_info[child] = i;
    for (lf::base::size_type n = 0; n < num_cells; ++n)
      is_part_of[lvl - 1][n] = parent_info[is_part_of[lvl][n]];
  }

  // Build new function objects operating on the finest mesh instead of one of
  // the coarser ones
  std::map<lf::base::size_type, std::function<Ret(const lf::mesh::Entity &,
                                                  const Eigen::Vector2d &)>>
      f_fine;
  const auto mesh_finest = meshes.getMesh(meshes.NumLevels() - 1);
  for (const auto [lvl, func] : f) {
    const auto mesh = meshes.getMesh(lvl);
    const auto lvl_is_part_of = is_part_of[lvl];
    auto fine_func = [mesh, mesh_finest, lvl_is_part_of, func = func](
                         const lf::mesh::Entity &entity,
                         const Eigen::Vector2d &x) -> Ret {
      const lf::base::size_type idx_finest = mesh_finest->Index(entity);
      const lf::base::size_type idx_coarse = lvl_is_part_of[idx_finest];
      const auto entity_coarse = mesh->EntityByIndex(0, idx_coarse);
      return func(*entity_coarse, x);
    };
    f_fine[lvl] = fine_func;
  }

  return f_fine;
}

}  // end namespace post_processing

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_POST_PROCESSING_MESH_HIERARCHY_FUNCTION_H
