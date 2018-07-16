#ifndef _LF_REFINEMENT_HIER_H_
#define _LF_REFINEMENT_HIER_H_

/**
 * @file refinement_hierarchy.h
 * @brief Data structures for managing nested meshes created by refinement
 *
 */

#include "refinement.h"

namespace lf::refinement {

/**
 * @brief Information about the refinement status of a point
 *
 * This information is used by methods of the class `MeshHierarchy` to
 * organize uniform and adaptive refinement.
 *
 * The key pieces of information are
 * - a flag indicating whether the point is to be duplicated (rp_copy)
 * - global index of the child point
 */
  struct PointChildInfo {
    PointChildInfo(void):ref_pat_(RefPat::rp_nil),child_point_idx_(-1) {}
    RefPat ref_pat_;
    lf::base::glb_idx_t child_point_idx_;
  };

/**
 * @brief Information about the refinement status of an edge
 *
 * This information is used by methods of the class `MeshHierarchy` to
 * organize uniform and adaptive refinement.
 *
 * The key pieces of information are
 * - the topological refinement type (rp_copy or rp_split)
 * - indices of all _interior_ child entities (edges or points)
 */
  struct EdgeChildInfo {
    EdgeChildInfo(void):ref_pat_(RefPat::rp_nil) {}
    RefPat ref_pat_;
    std::vector<lf::base::glb_idx_t> child_edge_idx_;
    std::vector<lf::base::glb_idx_t> child_point_idx_;
  };

/**
 * @brief Information about the refinement status of a cell
 *
 * This information is used by methods of the class `MeshHierarchy` to
 * organize uniform and adaptive refinement.
 *
 * The key pieces of information are
 * - the topological refinement type along with the anchor index used
 *   in the case of a cell, which completely defines the geometric refinement.
 * - indices of all _interior_ child entities
 */
  struct CellChildInfo {
    CellChildInfo(void):ref_pat_(RefPat::rp_nil), anchor_(-1) {}
    RefPat ref_pat_;
    lf::base::sub_idx_t anchor_;
    std::vector<lf::base::glb_idx_t> child_cell_idx_;
    std::vector<lf::base::glb_idx_t> child_edge_idx_;
    std::vector<lf::base::glb_idx_t> child_point_idx_;
  };
  
/**
 * @brief Information about possible parent entities
 */
struct ParentInfo {
  explicit ParentInfo(void) : parent_ptr_(nullptr), child_number_(-1) {}
  // Data members
  const mesh::Entity
      *parent_ptr_; /**< parent entity, not necessarily the same type */
  lf::base::sub_idx_t child_number_; /**< local index in the parent entity */
};

/**
 * @brief A hierarchy of nested 2D hybrid meshes created by refinement
 */
class MeshHierarchy {
 public:
  /**
   * @brief Initialize mesh hierarchy with an existing coarsest mesh
   *
   * @param base_mesh valid pointer to coarsest mesh
   * @param mesh_factory factory object creating new meshes during refinement
   */
  MeshHierarchy(std::shared_ptr<mesh::Mesh> base_mesh,
                mesh::MeshFactory &mesh_factory);
  MeshHierarchy(const MeshHierarchy &) = delete;
  MeshHierarchy &operator=(const MeshHierarchy &) = delete;

  /**
   * @brief number of meshes contained in the hierarchy, 1 for a single mesh
   */
  lf::base::size_type numLevels(void) const { return meshes_.size(); }

  /**
   * @brief access the mesh on a particular level
   */
  const mesh::Mesh &getMesh(lf::base::size_type level) {
    LF_VERIFY_MSG(level < meshes_.size(),
                  "Level " << level << " outside scope");
    return *meshes_.at(level);
  }

  /**
   * @brief Perform regular uniform refinement of the finest mesh in the
   * hierarchy
   *
   * This method carries out regular refinement of all cells of a mesh according
   * to the `rp_regular` refinement pattern.
   *
   * Regular refinement means that every node is copied, every edge is split
   * and every cell is subdivided into four smaller ones of the same shape.
   */
  void RefineRegular(void);
  /**
   * @brief Mark the edges of a mesh based on a predicate
   *
   * @param marker this should be functor of type
   * `std::function<bool(const Mesh &,const Entity &)>`
   * returning true if the passed edge is to be marked.
   */
  template <typename Marker>
  void MarkEdges(Marker &&marker);
  /**
   * @brief Conduct local refinement of the mesh splitting all marked edges
   */
  void RefineMarked(void);
  /**
   * @brief Destroy the mesh on the finest level unless it is the base mesh
   */
  void Coarsen(void);

  virtual ~MeshHierarchy(void) = default;

 private:
  /**
   * @brief Create new mesh according to refinement pattern
   *        provided for entities
   *
   * This method assumes that a refinement pattern has already be set
   * in the `ChildInfo` structure for each entity. According to this
   * information, refinement is carried out.
   */
  void PerformRefinement(void);

 private:
  /** the meshes managed by the MeshHierarchy object */
  std::vector<std::shared_ptr<mesh::Mesh>> meshes_;
  /** The mesh factory to be used to creating a new mesh */
  mesh::MeshFactory &mesh_factory_;
  /** information about children for each level and each class of entities */
  std::vector<std::vector<PointChildInfo>> point_child_infos_;
  std::vector<std::vector<EdgeChildInfo>> edge_child_infos_;
  std::vector<std::vector<CellChildInfo>> cell_child_infos_;
  /** information about parent entities on each level */
  std::vector<std::array<std::vector<ParentInfo>, 3>> parent_infos_;
  /** Information about marked edges */
  std::vector<std::vector<bool>> edge_marked_;
  /** Information about local refinement edges of triangles */
  std::vector<std::vector<lf::base::sub_idx_t>> refinement_edges_;
};

}  // namespace lf::refinement

#endif
