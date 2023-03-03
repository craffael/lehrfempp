#ifndef INCG_LF_REFINEMENT_HIER_H_
#define INCG_LF_REFINEMENT_HIER_H_

/**
 * @file
 * @brief Data structures for managing nested meshes created by refinement
 *
 */

#include <iostream>

#include "hybrid2d_refinement_pattern.h"

namespace lf::refinement {

/**
 * @brief Information about the refinement status of a point
 *
 * This information is used by methods of the class `MeshHierarchy` to
 * organize uniform and adaptive refinement.
 *
 * For a point there can be only two refinement states: non-refined (rp_nil) of
 * refined (rp_copy).
 *
 */
struct PointChildInfo {
  explicit PointChildInfo() = default;
  /** @brief a flag indicating whether the point is to be duplicated (rp_copy)
   */
  RefPat ref_pat{RefPat::rp_nil};
  /** @brief global index of the child point */
  glb_idx_t child_point_idx{idx_nil};
};

/**
 * @brief Information about the refinement status of an edge
 *
 * This information is used by methods of the class `MeshHierarchy` to
 * organize uniform and adaptive refinement.
 *
 * The key pieces of information are
 * - the topological refinement type (rp_copy or rp_split beside rp_nil in the
 *   case of non refinement)
 * - indices of all _interior_ child entities (edges or points). Their numbering
 *   advaces from enpoint #0 to endpoint #1.
 *
 * @sa MeshHierarchy::PerformRefinement()
 */
struct EdgeChildInfo {
  explicit EdgeChildInfo() = default;
  /** @brief type of refinement edge has undergone, see `RefPat` */
  RefPat ref_pat_{RefPat::rp_nil};
  /** @brief global indices of child edges in fine mesh */
  std::vector<glb_idx_t> child_edge_idx;
  /** @brief global indices of _interior_ child points in fine mesh */
  std::vector<glb_idx_t> child_point_idx;
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
 *
 * @sa MeshHierarchy::PerformRefinement()
 */
struct CellChildInfo {
  explicit CellChildInfo() = default;
  RefPat ref_pat_{RefPat::rp_nil};
  sub_idx_t anchor_{idx_nil};
  std::vector<glb_idx_t> child_cell_idx;
  std::vector<glb_idx_t> child_edge_idx;
  std::vector<glb_idx_t> child_point_idx;
};

/**
 * @brief Information about possible parent entities
 *
 * Also contains "relative geometry information" with respect to the
 * reference coordinate system of the parent entity.
 */
struct ParentInfo {
  explicit ParentInfo() = default;
  // Data members
  const mesh::Entity *parent_ptr{
      nullptr}; /**< parent entity, not necessarily the same type */
  glb_idx_t parent_index{
      idx_nil}; /**< index of parent entity w.r.t. coarse mesh */
  sub_idx_t child_number{idx_nil}; /**< local index in the parent entity */
  /** geometry object describing shape in parent entity in relative coordinates
   */
  std::unique_ptr<const lf::geometry::Geometry> rel_ref_geo_{nullptr};
};

/**
 * @brief A hierarchy of nested 2D hybrid meshes created by refinement
 *
 * This is the fundamental class for managing sequences of meshes that have been
 * created by successive local or global refinement of an initial coarsest mesh.
 * Objects of this class store comprehensive information about the relationship
 * of the meshes in the sequence: parent-child relationships of mesh cells and
 * local refinement patterns.
 *
 * #### sample usage
 * @snippet meshhierarchy.cc usage
 *
 * In the context of studying convergence of finite elemment method objects of
 * type MeshHierarchy can be employed to hold sequencies of uniformly refine
 * meshes, see @ref MeshHierarchy::RefineRegular().
 */
class MeshHierarchy {
 public:
  /**
   * @brief Initialize mesh hierarchy with an existing coarsest mesh
   *
   * @param base_mesh valid pointer to _non-const_ coarsest mesh
   * @param mesh_factory factory object creating new meshes during refinement
   *
   * - Stores shared pointer to coarsest mesh.
   * - Sets _refinement edges_ of all cells according to the
   * longest-edge criterion.
   * - Initializes `ChildInfo` data structures of all entities to indicate
   * absence of children, since no refinement has been done yet.
   * - Fills void `ParentInfo` data structure, since the mesh has not been
   * created by refinement.
   * - Unmarks all edges, see RefineMarked().
   */
  MeshHierarchy(const std::shared_ptr<mesh::Mesh> &base_mesh,
                std::unique_ptr<mesh::MeshFactory> mesh_factory);
  MeshHierarchy(const MeshHierarchy &) = delete;
  MeshHierarchy &operator=(const MeshHierarchy &) = delete;
  MeshHierarchy(MeshHierarchy &&) = delete;
  MeshHierarchy &operator=(MeshHierarchy &&) = delete;

  /**
   * @brief number of meshes contained in the hierarchy, 1 for a single mesh
   */
  [[nodiscard]] size_type NumLevels() const { return meshes_.size(); }

  /**
   * @brief access the mesh on a particular level
   *
   * @param level specifies level of interest, 0 stands for coarsest level
   * @return shared pointer to mesh on specified level
   */
  [[nodiscard]] std::shared_ptr<const mesh::Mesh> getMesh(
      size_type level) const {
    LF_VERIFY_MSG(level < meshes_.size(),
                  "Level " << level << " outside scope");
    return meshes_.at(level);
  }
  /**
   * @copydoc lf::refinement::MeshHierarchy::getMesh()
   */
  [[nodiscard]] std::shared_ptr<mesh::Mesh> getMesh(size_type level) {
    LF_VERIFY_MSG(level < meshes_.size(),
                  "Level " << level << " outside scope");
    return meshes_.at(level);
  }

  /**
   * @brief Provides array of shared pointers to meshes contained in the
   * hierarchy
   *
   * @return vector of shared pointers to lf::mesh:Mesh objects
   */
  [[nodiscard]] std::vector<std::shared_ptr<const mesh::Mesh>> getMeshes()
      const {
    return {meshes_.begin(), meshes_.end()};
  }

  /**
   * @brief Obtain refinement information for all points
   *
   * @param level refinement level to be queried
   * @return vector for PointChildInfo record for every node
   *
   * @sa PointChildInfo
   */
  [[nodiscard]] const std::vector<PointChildInfo> &PointChildInfos(
      size_type level) const {
    LF_VERIFY_MSG(level < NumLevels(), "Illegal level " << level);
    return point_child_infos_[level];
  }
  /**
   * @brief Obtain refinement information for all edges
   *
   * @param level refinement level to be queried
   * @return vector for EdgeChildInfo record for every node
   *
   * @sa EdgeChildInfo
   */
  [[nodiscard]] const std::vector<EdgeChildInfo> &EdgeChildInfos(
      size_type level) const {
    LF_VERIFY_MSG(level < NumLevels(), "Illegal level " << level);
    return edge_child_infos_[level];
  }
  /**
   * @brief Obtain refinement information for all
   *
   * @param level refinement level to be queried
   * @return vector for CellChildInfo record for every node
   *
   * @sa CellChildInfo
   */
  [[nodiscard]] const std::vector<CellChildInfo> &CellChildInfos(
      size_type level) const {
    LF_VERIFY_MSG(level < NumLevels(), "Illegal level " << level);
    return cell_child_infos_[level];
  }
  /**
   * @brief Fetch information about parents
   *
   * @param level refinement level of interest
   * @param codim co-dimension of entities to be queried
   * @return vector for ParentInfo record for every entity
   *
   * @sa ParentInfo
   */
  [[nodiscard]] const std::vector<ParentInfo> &ParentInfos(size_type level,
                                                           dim_t codim) const {
    LF_VERIFY_MSG(level < NumLevels(), "Illegal level " << level);
    LF_VERIFY_MSG(codim < 3, "Codim = " << codim << " illegal");
    return parent_infos_[level][codim];
  }
  /**
   * @brief Access refinement edge indices
   *
   * @param level refinement level of interest
   * @return vector of (local) sub-entity index of refinement edge for every
   * cell
   */
  [[nodiscard]] const std::vector<sub_idx_t> &RefinementEdges(
      size_type level) const {
    LF_VERIFY_MSG(level < NumLevels(), "Illegal level " << level);
    return refinement_edges_[level];
  }

  /**
   * @brief Perform regular or barycentric uniform refinement of the finest mesh
   * in the hierarchy
   *
   * @param ref_pat selector for type of uniform refinement: default is
   * rp_regular, rp_barycentric choses barycentric refinement.
   *
   * This method carries out uniform refinement of all cells of a mesh according
   * to the `rp_regular` or `rp_barycentric` refinement patterns.
   *
   * A new mesh is added to the bottom of the hierarchy by regularly refining
   * the _finest mesh_ in the hierarchy.
   * Regular refinement means that every node is copied, every edge is split
   * and every cell is subdivided into four or six smaller ones of the same
   * shape.
   *
   * Internally, this method flags all nodes as to be copied, all edges as to be
   * split and all cells as to be refined according to the passed refinement
   * pattern. Then it calls PerformRefinement().
   */
  void RefineRegular(RefPat ref_pat = RefPat::rp_regular);
  /**
   * @brief Mark the edges of a mesh based on a predicate
   *
   * @param marker this should be functor of type
   * `std::function<bool(const Mesh &,const Entity &)>`
   * returning true if the passed edge is to be marked.
   *
   * The _marker_ object also takes a reference to a mesh, because
   * marking makes sense only for the finest level. The mesh on the
   * finest level is provided to the marker object by the `MeshHierarchy`.
   *
   * Of course, marking will always affect the finest mesh in hierarchy.
   */
  template <typename Marker>
  void MarkEdges(Marker &&marker);

  /**
   * @brief Conduct local refinement of the mesh splitting all marked edges
   *
   * This method creates a new mesh by selectively (locally) refining entities
   * of the current finest mesh in the hierarchy. Refinement is controlled by
   * the boolean vector edge_marked_  that indicates, which edges must be
   * refined (= split) in the course of refinement.
   *
   * ### Algorithm
   *
   * - First all marked edges are labelled as "to be split".
   * - REPEAT
   *   + Set refinement pattern of all cells to accommodate edges to be split
   *   + Add "to be split" tag to edges according to local refinement pattern
   * for cells
   *
   *   UNTIL no _extra_ edges had to be tagged as "to be split"
   *
   * For details please consult the comments in mesh_hierarchy.cc
   *
   * This algorithm ends with a set of local refinement patterns for every
   * entity that is compatible with a _conforming_ finite element mesh, that is,
   * hanging nodes are avoided.
   */
  void RefineMarked();
  /**
   * @brief _Destroy_ the mesh on the finest level unless it is the base mesh
   *
   * @note the use of shared pointers prevents destruction if the finest mesh
   *       is still in use somewhere else in the code.
   */
  void Coarsen();
  /**
   * @brief shape of child entity in parent's reference coordinates
   *
   * @param level level of the fine mesh in the mesh hierarchy
   * @param e reference to entity object of the fine mesh
   * @return pointer to geometry description of entity of fine mesh with respect
   * to reference coordinate system of its parent entity.
   *
   * We consider the shape of the reference entity of the parent entity of e as
   * a triangulated domain, if the parent entity has a positive dimension, that
   * is, if it is not a point. The entities of that triangulation correspond to
   * the reference shapes of the child entities on the fine mesh. This method
   * returns the shape of an entity of that triangulation.
   *
   * Example: The case of a quadrilateral split into four child entities,
   * corresponding to `rp_threeedge` as defined in @ref lf::refinement::RefPat,
   * @ref lf::refinement::Hybrid2DRefinementPattern::ChildPolygons().
   * The following image gives the local child numbers in pink.
   * @image html refinement_quad/rp_threeedge_unitsquare.png width=500px
   * The geometries of the child entities are as follows:
   *
   * - Child 0: quadrilateral with corners
   * \f[ \left[\begin{array}{cccc} 1 & 0 & 0 & 1 \\ 1 & 1 & 0.5 & 0.5
   * \end{array}\right] \f]
   * - Child 1: triangle with corners
   * \f[ \left[\begin{array}{ccc} 0 & 0.5 & 0 \\ 0 & 0  & 0.5 \end{array}
   * \right] \f]
   * - Child 2: triangle with corners
   * \f[ \left[\begin{array}{ccc} 1 & 0.5 & 1 \\ 0.5 & 0 & 0 \end{array} \right]
   * \f]
   * - Child 3: triangle with corners
   * \f[ \left[\begin{array}{ccc} 0.5 & 1 & 0 \\ 0 & 0.5 & 0.5 \end{array}
   * \right] \f]
   *
   * One of the  displayed matrices is returned if the entity `e` corresponds to
   * that child of a quadrilateral.
   */
  [[nodiscard]] const lf::geometry::Geometry *GeometryInParent(
      size_type level, const lf::mesh::Entity &e) const;

  /** @brief Retrieve the parent of an entity contained in a mesh of a
   * refinement hierarchy
   *  @param level refinement level > 0 of the fine mesh
   *  @param e pointer to (potential) child entity, must be contained in mesh on
   * the specified level.
   *  @return pointer to parent entity
   *
   * @note to be used in combination with geometry::Geometry::GeometryInParent()
   */
  [[nodiscard]] const lf::mesh::Entity *ParentEntity(
      size_type level, const lf::mesh::Entity &e) const;

  /**
   * @brief Output of information about the mesh hierarchy.
   *
   * @param o output stream, can be `std::cout` or similar
   * @param ctrl controls the level of detail of the generated output.
   * @return output stream
   *
   * The type of output is controlled by the `ctrl` static control
   * variable. If its second bit is set, the output function of the
   * mesh class is used.
   *
   * This is a rudimentary implementation and should be extended.
   *
   */
  std::ostream &PrintInfo(std::ostream &o, unsigned ctrl = 0) const;

  virtual ~MeshHierarchy() = default;

 private:
  /**
   * @brief Create new mesh according to refinement pattern
   *        provided for entities
   *
   * This function expects that the refinement patterns  stored in the
   * vectors `point_child_infos_`, `edge_child_infos_` and `cell_child_infos_`
   * have been initialized consistently for the finest mesh. According to this
   * information, refinement is carried out using the object pointed to by
   * mesh_factory_ to created new entities by calling
   * @ref lf::mesh::MeshFactory::Build().
   *
   * The vectors  `point_child_infos_`, `edge_child_infos_` and
   * `cell_child_infos_` will be augmented with information about the indices of
   * the child entities contained  in the newly created finest mesh.
   *
   * This method relies on @ref lf::geometry::Geometry::ChildGeometry() to
   * obtain information about the shape of child entities in the form of
   * lf::geometry::Geometry objects.
   *
   * The method also initializes the data vectors in `_parent_infos_` for
   * the newly created finest mesh.
   *
   */
  void PerformRefinement();

  /**
   * @brief Initialization of rel_ref_geo fields of @ref ParentInfo structures
   *
   * This method created @ref geometry::Geometry type objects describing
   * the shape of a child entity in the reference coordinates of its parent.
   *
   * @note This method assumes that the parent-child connections of the mesh
   * hierarchy have been initialized completely already. Therefore this method
   * is invoked at the end of @ref PerformRefinement().
   */
  void initGeometryInParent();

  /** @brief the meshes managed by the MeshHierarchy object */
  std::vector<std::shared_ptr<mesh::Mesh>> meshes_;
  /** @brief The mesh factory to be used to creating a new mesh */
  std::unique_ptr<mesh::MeshFactory> mesh_factory_;
  /** @brief information about children of nodes for each level */
  std::vector<std::vector<PointChildInfo>> point_child_infos_;
  /** @brief information about children of edges for each level */
  std::vector<std::vector<EdgeChildInfo>> edge_child_infos_;
  /** @brief information about children of cells for each level */
  std::vector<std::vector<CellChildInfo>> cell_child_infos_;
  /** @brief information about parent entities on each level */
  std::vector<std::array<std::vector<ParentInfo>, 3>> parent_infos_;
  /** @brief Information about marked edges */
  std::vector<std::vector<bool>> edge_marked_;
  /** @brief Information about local refinement edges of triangles */
  std::vector<std::vector<sub_idx_t>> refinement_edges_;

  /**
   * @brief Finds the index of the longest edge of a triangle
   *
   * This method is used for setting refinement edges on coarsest meshes.
   * Called in the constructor of MeshHierarchy.
   */
  [[nodiscard]] static sub_idx_t LongestEdge(const lf::mesh::Entity &T);

 public:
  /**
   * @brief Is used by MeshHierarchy to log additional information for debugging
   * purposes.
   */
  static std::shared_ptr<spdlog::logger> &Logger();
};

template <typename Marker>
void MeshHierarchy::MarkEdges(Marker &&marker) {
  // Retrieve the finest mesh in the hierarchy
  const mesh::Mesh &finest_mesh(*meshes_.back());

  LF_VERIFY_MSG(edge_marked_.back().size() == finest_mesh.NumEntities(1),
                "Length  mismatch for edge flag array");

  // Run through the edges = entities of co-dimension 1
  for (const mesh::Entity *edge : finest_mesh.Entities(1)) {
    glb_idx_t edge_index = finest_mesh.Index(*edge);
    (edge_marked_.back())[edge_index] = marker(finest_mesh, *edge);
  }
}  // end MeshHierarchy::MarkEdges

/**
 * @brief Generated a sequence of nested 2D hybrid mehes by regular or
 * barycentric refinement
 *
 * @param mesh_p pointer to the coarsest mesh, from which refinment will start
 * @param ref_lev desired number of refinement steps
 * @param ref_pat at selector for type of uniform refinement: default is
 * rp_regular, rp_barycentric choses barycentric refinement.
 * @return shared pointer to a MeshHierarchy object.
 *
 * Relies on lf::mesh::hybrid2d::MeshFactory as builder class for mesh entities.
 * Invokes the method MeshHierarhy::RefineRegular() for refinement.
 */
std::shared_ptr<MeshHierarchy> GenerateMeshHierarchyByUniformRefinemnt(
    const std::shared_ptr<lf::mesh::Mesh> &mesh_p, lf::base::size_type ref_lev,
    RefPat ref_pat = RefPat::rp_regular);

/**
 * @brief Utility class: selection of entities according to the position of
 * their midpoint
 * @tparam POSPRED predicate depending on physical point location
 *
 * ### Type requirements
 * POSPRED must have an evaluation operator that accepts an Eigen matrix/vector
 * object whose columns are viewed as coordinate vectors. It should return a
 * boolean value.
 */
template <typename POSPRED>
class EntityCenterPositionSelector {
 public:
  EntityCenterPositionSelector(const EntityCenterPositionSelector &) = default;
  EntityCenterPositionSelector(EntityCenterPositionSelector &&) noexcept =
      default;
  EntityCenterPositionSelector &operator=(
      const EntityCenterPositionSelector &) = default;
  EntityCenterPositionSelector &operator=(
      EntityCenterPositionSelector &&) noexcept = default;
  /**
   * @brief Preparing reference coordinates of "centers"
   *
   * @param pos_pred object for true/false classification of physicals points
   */
  explicit EntityCenterPositionSelector(POSPRED pos_pred)
      : pos_pred_(std::move(pos_pred)) {}
  /** @brief Operator testing location of "center"
   *  @param ent reference to a mesh entity
   *  @return value of the predicate when given the center of the entity
   */
  bool operator()(const lf::mesh::Entity &ent) const {
    const lf::base::RefEl ref_el_type = ent.RefEl();
    // Obtain shape of entity
    const lf::geometry::Geometry *geo_ptr = ent.Geometry();
    LF_ASSERT_MSG(geo_ptr != nullptr, "Missing geometry for " << ent);
    switch (ref_el_type) {
      case lf::base::RefEl::kPoint(): {
        const Eigen::MatrixXd pos(geo_ptr->Global(kpoint_center_));
        return pos_pred_(pos);
      }
      case lf::base::RefEl::kSegment(): {
        const Eigen::MatrixXd pos(geo_ptr->Global(kedge_center_));
        return pos_pred_(pos);
      }
      case lf::base::RefEl::kTria(): {
        const Eigen::MatrixXd pos(geo_ptr->Global(ktria_center_));
        return pos_pred_(pos);
      }
      case lf::base::RefEl::kQuad(): {
        const Eigen::MatrixXd pos(geo_ptr->Global(kquad_center_));
        return pos_pred_(pos);
      }
      default: {
        LF_ASSERT_MSG(false, "Illegal entity type");
        break;
      }
    }  // end switch
    return false;
  }

  virtual ~EntityCenterPositionSelector() = default;

 private:
  /** object for true/false classification of of physicals points */
  POSPRED pos_pred_;

  static const Eigen::MatrixXd kpoint_center_;
  static const Eigen::MatrixXd kedge_center_;
  static const Eigen::MatrixXd ktria_center_;
  static const Eigen::MatrixXd kquad_center_;
};

template <typename POSPRED>
const Eigen::MatrixXd EntityCenterPositionSelector<POSPRED>::kpoint_center_ =
    Eigen::MatrixXd::Zero(0, 1);

template <typename POSPRED>
const Eigen::MatrixXd EntityCenterPositionSelector<POSPRED>::kedge_center_ =
    (Eigen::MatrixXd(1, 1) << 0.5).finished();

template <typename POSPRED>
const Eigen::MatrixXd EntityCenterPositionSelector<POSPRED>::ktria_center_ =
    (Eigen::MatrixXd(2, 1) << 0.33, 0.33).finished();

template <typename POSPRED>
const Eigen::MatrixXd EntityCenterPositionSelector<POSPRED>::kquad_center_ =
    (Eigen::MatrixXd(2, 1) << 0.5, 0.5).finished();

}  // namespace lf::refinement

#endif
