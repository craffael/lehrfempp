#ifndef _LF_DOFHD_H
#define _LF_DOFHD_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief A general DOF handler interface
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include <lf/mesh/mesh.h>
#include "assembly_types.h"

namespace lf::assemble {
/**
 * @brief A general (_interface_) class for DOF handling.
 *
 * Objects of this class provide the local-to-global map for indices
 * of local/global shape functions.
 *
 * # Rules for ordering global shape functions <-> degrees of freedom (dof)
 *
 * Note that every global shape (dof) function belongs to a unique mesh entity.
 * Ordering them means assigning a unique (vector component) index starting from
 * zero until `NoDofs-1`.
 * -# Dofs belonging to entities of a larger co-dimension are arranged before
 *    dofs associated with entities of a smaller co-dimension:
 *
 *  This means **First dofs on nodes, then dofs on edges, then dofs on cells**
 * -# Dofs owned by the same mesh entity are numbered contiguously.
 *
 * -# if two entities have the same co-dimension, then the dofs of that
 *    with the lower index are numbered before the dofs of that with
 *    the higher index: *dof numbering is compatible with indexing*
 *
 * # Rules for numbering local shape functions (local dof)
 *
 * -# As above, Dofs belonging to _sub-entities_ of a larger (relative)
 * co-dimension are arranged before dofs associated with _sub-entities_ of a
 * smaller co-dimension.
 * -# As above, dofs belonging to the same sub-entity are numbered contiguously.
 * -# dofs for _sub-entities_ of the same co-dimension are taken into account in
 * the order given by the _local indexing_ of the sub-entities.
 *
 * # Local-to-global dof index mapping for entities of higher co-dimension
 *
 * The method getGlobalDofs() is available for entities of _any co-dimension_.
 * For instance, this allows assembly of contributions from low-dimensional
 * manifolds like boundaries or interfaces.
 *
 * However, the consistency of the local numbering of dofs has to be ensured.
 * In two dimensions this is an issue for edges only and can be resolved by
 * taking into account the orientation (direction) of the edge: dofs are ordered
 * along the edge and those "closer to endpoint 0" are numbered first. In 3D
 * many more situations have to be dealt with.
 */
class DofHandler {
 public:
  /**
   * @brief The constructor just stores the mesh pointer
   */
  DofHandler() = default;
  virtual ~DofHandler() = default;

  /** Disabled constructors and assignment operators */
  /**@{*/
  /** @name Disabled */
  DofHandler(const DofHandler &) = delete;
  DofHandler(DofHandler &&) = delete;
  DofHandler &operator=(const DofHandler &) = delete;
  DofHandler &operator=(DofHandler &&) = delete;
  /**@}*/

  /**
   * @brief total number of dof's handled by the object
   */
  virtual size_type NoDofs() const = 0;

  /**
   * @brief tells the number of degrees of freedom subordinated to an entity
   *
   * @param entity reference to an entity of the underlying mesh
   *
   * This method informs about the length of the vector returned by
   * GlobalDofIndices().
   * @sa GlobalDofIndices()
   */
  virtual size_type NoLocalDofs(const lf::mesh::Entity &entity) const = 0;

  /**
   * @brief provides number of shape functions associated with an entity
   *
   * @param entity entity of underlying mesh whose number of shape functions
   *               is queried
   *
   * The return value of this method must be equal to the length of the range
   * built by GlobalDofIndices().
   */
  virtual size_type NoInteriorDofs(const lf::mesh::Entity &entity) const = 0;

  /**
   * @brief access to indices of global dof's belonging to an entity
   *
   * @param entity reference to the entity for which the dof's are to be
   *        fetched. This entity must belong to the underlying mesh.
   * @return cardinal number range of global dof indices
   *
   * The basis functions of every finite element space must be associated with a
   * unique geometric entity. Conversely, every entity can possess a finite
   * number of basis functions = degrees of freedom. This functions return the
   * global indices of all basis functions associated with the entity _and its
   * sub-entitites_.
   *
   * The size of the returned range must agree with the value returned
   * by NoLocalDofs() when supplied with the same arguments.
   */
  virtual lf::base::RandomAccessRange<const gdof_idx_t> GlobalDofIndices(
      const lf::mesh::Entity &entity) const = 0;

  /**
   * @brief global indices of shape functions _belonging to an entity_
   *
   * @param entity entity for which shape functin indices are queried
   * @return cardinal number range of global indices of shape functions
   *
   * Each global shape function is associated with a unique mesh entity.
   * This method provides all the global indices of the shape function
   * associated to the entity specified by the function arguments.
   */
  virtual lf::base::RandomAccessRange<const gdof_idx_t>
  InteriorGlobalDofIndices(const lf::mesh::Entity &entity) const = 0;

  /**
   * @brief retrieve unique entity at which a basis function is located
   *
   * @param dofnum global index of the basis function/degree of freedom
   *
   * This function returns the unique geometric entity to which a particular
   * basis function = degree of freedom is associated.
   *
   * This function is hardly ever needed in finite element codes and is supplied
   * for debugging purposes.
   * @sa GlobalDofIndices()
   */
  virtual const lf::mesh::Entity &Entity(gdof_idx_t dofnum) const = 0;

  /** @brief Acess to underlying mesh object
   *
   * Every DofHandler object has to be associated with a unique mesh.
   * All entities passed to the DofHandler object must belong to that mesh.
   */
  virtual std::shared_ptr<const lf::mesh::Mesh> Mesh() const = 0;

  /** @brief output control variable
   *
   * This variable is used by the output operator for DofHandler objects.
   *
   * Initialized with 0
   *
   * The following output levels are supported
   * - if = 0: minimal out: just prints global no of dofs
   * - if > 0 and divisible by 2: print dofs associated with all entities
   * - if > 0 and divisible by 3: print entities belonging to dofs
   * - if > 0 and divisible by 10: also print interior dof indices
   */
  static int output_ctrl_;
};

/** @brief output operator for DofHandler objects */
std::ostream &operator<<(std::ostream &o, const DofHandler &dof_handler);

/* ====================================================================== */

/**
 * @brief Dofhandler for uniform finite element spaces
 *
 * This management class for indices of global shape functions
 * is suitable for situations where every geometric entity of a particular
 * type has exactly the same number of shape functions belonging to it.
 */
class UniformFEDofHandler : public DofHandler {
 public:
  /** Construction from local dof layout */
  /**@{*/
  /** @name Constructors */

  /** @brief Construction from a map object
   *
   * @param dofmap map telling number of interior dofs for every type of entity
   */
  using dof_map_t = std::map<lf::base::RefEl, base::size_type>;
  UniformFEDofHandler(std::shared_ptr<const lf::mesh::Mesh> mesh,
                      dof_map_t dofmap);
  /**@}*/

  /**
   * @copydoc DofHandler::GetNoDofs()
   */
  size_type NoDofs() const override { return num_dof_; }

  /**
   * @copydoc DofHandler::GetNoLocalDofs()
   * @sa GlobalDofIndices()
   */
  size_type NoLocalDofs(const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GetNoInteriorDofs()
   * @sa InteriorGlobalDofIndices
   */
  size_type NoInteriorDofs(const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GlobalDofIndices()
   */
  lf::base::RandomAccessRange<const gdof_idx_t> GlobalDofIndices(
      const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::InteriorGlobalDofIndices()
   */
  /** @copydoc DofHandler::InteriorGlobalDofIndices() */
  lf::base::RandomAccessRange<const gdof_idx_t> InteriorGlobalDofIndices(
      const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GetEntity()
   * @sa GlobalDofIndices()
   */
  const lf::mesh::Entity &Entity(gdof_idx_t dofnum) const override {
    LF_VERIFY_MSG(dofnum < dof_entities_.size(),
                  "Illegal dof index " << dofnum << ", max = " << num_dof_);
    return *dof_entities_[dofnum];
  }

  /** @copydoc DofHandler::mesh()
   */
  std::shared_ptr<const lf::mesh::Mesh> Mesh() const override { return mesh_; }

 private:
  /**
   * @brief initialization of internal index arrays
   */
  void initIndexArrays();
  /** @brief compute number of shape functions covering an entity type
   *
   * This method assumes that the variables  no_loc_dof_point_,
   * no_loc_dof_segment_, no_loc_dof_tria_, no_loc_dof_quad_
   * have already been set.
   *
   * @sa LocalStaticDOFs2D::TotalNoLocDofs()
   */
  void initTotalNoDofs();

  // Access method to numbers and values of indices of shape functions
  lf::base::RandomAccessRange<const gdof_idx_t> GlobalDofIndices(
      lf::base::RefEl ref_el_type, glb_idx_t entity_index) const;

  lf::base::RandomAccessRange<const gdof_idx_t> InteriorGlobalDofIndices(
      lf::base::RefEl ref_el_type, glb_idx_t entity_index) const;

  size_type GetNoLocalDofs(lf::base::RefEl ref_el_type,
                           glb_idx_t /*unused*/) const {
    return NoCoveredDofs(ref_el_type);
  }

  /** co-dimensions of different geometric entities in a 2D mesh */
  /**@{*/
  const size_type kNodeOrd = 2; /**< */
  const size_type kEdgeOrd = 1; /**< */
  const size_type kCellOrd = 0; /**< */
  /**@}*/

  /** Number of covered dofs for an entity type */
  inline size_type NoCoveredDofs(lf::base::RefEl ref_el_type) const {
    size_type no_covered_dofs;
    switch (ref_el_type) {
      case lf::base::RefEl::kPoint(): {
        no_covered_dofs = no_dofs_[kNodeOrd];
        break;
      }
      case lf::base::RefEl::kSegment(): {
        no_covered_dofs = no_dofs_[kEdgeOrd];
        break;
      }
      case lf::base::RefEl::kTria(): {
        no_covered_dofs = num_dofs_tria_;
        break;
      }
      case lf::base::RefEl::kQuad(): {
        no_covered_dofs = num_dofs_quad_;
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Illegal entity type");
        break;
      }
    }  // end switch
    return no_covered_dofs;
  }

  /** Number of interior shape functions for an entity type */
  inline size_type NoInteriorDofs(lf::base::RefEl ref_el_type) const {
    size_type no_loc_dofs;
    switch (ref_el_type) {
      case lf::base::RefEl::kPoint(): {
        no_loc_dofs = no_loc_dof_point_;
        break;
      }
      case lf::base::RefEl::kSegment(): {
        no_loc_dofs = no_loc_dof_segment_;
        break;
      }
      case lf::base::RefEl::kTria(): {
        no_loc_dofs = no_loc_dof_tria_;
        break;
      }
      case lf::base::RefEl::kQuad(): {
        no_loc_dofs = no_loc_dof_quad_;
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Illegal entity type");
        break;
      }
    }  // end switch
    return no_loc_dofs;
  }

  /** The mesh on which the degrees of freedom are defined */
  std::shared_ptr<const lf::mesh::Mesh> mesh_;
  /** The total number of degrees of freedom */
  size_type num_dof_{0};
  /** Vector of entities to which global basis functions are associated */
  std::vector<const lf::mesh::Entity *> dof_entities_;
  /** Vectors of global indices of dofs belonging to entities of different
      topological type */
  std::array<std::vector<gdof_idx_t>, 3> dofs_;
  /** Number of dofs covering entities of a particular type */
  std::array<size_type, 3> no_dofs_;
  /** (Maximum) number of shape functions covering entities
   * of a particular co-dimension  */
  size_type num_dofs_tria_{0}, num_dofs_quad_{0};
  /** Numbers of local _interior_ shape functions associated with different
   * types of entities */
  /**@{*/
  size_type no_loc_dof_point_{0};
  size_type no_loc_dof_segment_{0};
  size_type no_loc_dof_tria_{0};
  size_type no_loc_dof_quad_{0};
  /**@}*/
};

/* ====================================================================== */

/** @brief Dof handler allowing variable local dof layouts
 *
 * This dof handler can accommodate cases where entities of the mesh
 * have different numbers of local shape functions attached to them.
 * This is relevant, for instance, for hp-FEM.
 *
 */
class DynamicFEDofHandler : public DofHandler {
 public:
  /** @brief Set-up of dof handler
   *
   * @tparam LOCALDOFINFO type for object telling number of interior local shape
   * functions
   * @param mesh_p pointer to underlying mesh
   * @param locdof functor object telling number of _interior_ dofs for every
   * entity of the mesh.
   *
   * This constructor performs the initialization of all internal index arrays.
   *
   * ### Type requirements for LOCALDOFINFO
   *
   * This type must provide method
   *
   *      size_type operator (const lf::mesh::Entity &);
   *
   * that returns the number of local shape functions associated with an entity.
   * Note that this gives the number of so-called interior local shape functions
   * belonging to an entity and not the number of local shape functions whose supports
   * will include the entity. Also refer to the documentation of DynamicFEDofHandler::NoLocalDofs()
   * and DynamicFEDofHandler::NoInteriorDofs(). 
   */
  template <typename LOCALDOFINFO>
  DynamicFEDofHandler(std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                      LOCALDOFINFO &&locdof)
      : mesh_p_(std::move(mesh_p)) {
    LF_ASSERT_MSG((mesh_p_->DimMesh() == 2), "Can handle 2D meshes only");

    // Index counter for global shape functions = global dof
    gdof_idx_t dof_idx = 0;

    // Step I: Set indices for shape functions on nodes
    // Run through node indices (entities of co-dimension 2)
    const size_type no_nodes = mesh_p_->Size(2);
    no_int_dofs_[2].resize(no_nodes, 0);
    offsets_[2].resize(no_nodes + 1, 0);
    // Traverse nodes (co-dimension-2 entities) based on indices
    for (glb_idx_t node_idx = 0; node_idx < no_nodes; node_idx++) {
      // Obtain pointer to node entity
      const mesh::Entity *node_p{mesh_p_->EntityByIndex(2, node_idx)};
      LF_ASSERT_MSG(mesh_p_->Index(*node_p) == node_idx, "Node index mismatch");
      // Offset for indices of node in index vector
      glb_idx_t node_dof_offset = dof_idx;
      offsets_[2][node_idx] = node_dof_offset;
      // Request number of local shape functions associated with the node
      size_type no_int_dof_node = locdof(*node_p);
      no_int_dofs_[2][node_idx] = no_int_dof_node;

      // Store dof indices in array
      for (int j = 0; j < no_int_dof_node; j++) {
        dofs_[2].push_back(dof_idx);
        dof_entities_.push_back(node_p);  // Store entity for current dof
        dof_idx++;                        // Move on to next index
      }
    }
    // Set sentinel
    offsets_[2][no_nodes] = dof_idx;

    // Step II: Set indices for shape functions on edges (co-dimension = 1)
    const size_type no_edges = mesh_p_->Size(1);
    // Set length of edge-related index vectors
    no_int_dofs_[1].resize(no_edges, 0);
    offsets_[1].resize(no_edges + 1, 0);
    // points to the position of the current dof index
    size_type edge_dof_offset = 0;
    for (glb_idx_t edge_idx = 0; edge_idx < no_edges; edge_idx++) {
      // Obtain pointer to edge entity
      const mesh::Entity *edge{mesh_p_->EntityByIndex(1, edge_idx)};
      LF_ASSERT_MSG(mesh_p_->Index(*edge) == edge_idx, "Edge index mismatch");
      // Offset for indices of edge dof in index vector
      offsets_[1][edge_idx] = edge_dof_offset;
      size_type no_int_dof_edge = locdof(*edge);
      no_int_dofs_[1][edge_idx] = no_int_dof_edge;

      // Obtain indices for basis functions sitting at endpoints
      // Endpoints are mesh entities with co-dimension = 2
      for (const lf::mesh::Entity &endpoint : edge->SubEntities(1)) {
        const glb_idx_t ep_idx(mesh_p_->Index(endpoint));
        glb_idx_t ep_dof_offset = offsets_[2][ep_idx];
        size_type no_int_dofs_ep = no_int_dofs_[2][ep_idx];
        // Copy indices of shape functions from nodes to edge
        for (int j = 0; j < no_int_dofs_ep; j++) {
          dofs_[1].push_back(dofs_[2][ep_dof_offset + j]);
          edge_dof_offset++;
        }
      }
      // Set indices for interior edge degrees of freedom
      for (int j = 0; j < no_int_dof_edge; j++) {
        dofs_[1].push_back(dof_idx);
        edge_dof_offset++;
        dof_entities_.push_back(edge);
        dof_idx++;
      }
    }
    // Set sentinel
    offsets_[1][no_edges] = edge_dof_offset;

    // Step III: Set indices for shape functions on cells (co-dimension = 0)
    const size_type no_cells = mesh_p_->Size(0);
    // Set length of cell-related index vectors
    no_int_dofs_[0].resize(no_cells, 0);
    offsets_[0].resize(no_cells + 1, 0);
    // points to the position of the current dof index
    size_type cell_dof_offset = 0;
    for (glb_idx_t cell_idx = 0; cell_idx < no_cells; cell_idx++) {
      // Obtain pointer to current ell
      const mesh::Entity *cell{mesh_p_->EntityByIndex(0, cell_idx)};
      // Offset for indices of cell dof in index vector
      offsets_[0][cell_idx] = cell_dof_offset;
      size_type no_int_dof_cell = locdof(*cell);
      no_int_dofs_[0][cell_idx] = no_int_dof_cell;

      // Obtain indices for basis functions in vertices
      for (const lf::mesh::Entity &vertex : cell->SubEntities(2)) {
        const glb_idx_t vt_idx(mesh_p_->Index(vertex));
        glb_idx_t vt_dof_offset = offsets_[2][vt_idx];
        size_type no_int_dofs_vt = no_int_dofs_[2][vt_idx];
        // Copy indices of shape functions from nodes to cell
        for (int j = 0; j < no_int_dofs_vt; j++) {
          dofs_[0].push_back(dofs_[2][vt_dof_offset + j]);
          cell_dof_offset++;
        }
      }

      // Collect indices of interior shape functions of edges
      // Internal ordering will depend on the orientation of the edge
      lf::base::RandomAccessRange<const lf::mesh::Orientation>
          edge_orientations(cell->RelativeOrientations());
      lf::base::RandomAccessRange<const lf::mesh::Entity> edges(
          cell->SubEntities(1));
      // Loop over edges
      for (int ed_sub_idx = 0; ed_sub_idx < cell->RefEl().NumSubEntities(1);
           ed_sub_idx++) {
        const glb_idx_t edge_idx = mesh_p_->Index(edges[ed_sub_idx]);
        const size_type no_int_dof_edge = no_int_dofs_[1][edge_idx];
        const glb_idx_t edge_int_dof_offset =
            offsets_[1][edge_idx + 1] - no_int_dof_edge;

        // Copy indices of shape functions from edges to cell
        // The order, in which they are copied depends on the relative
        // orientation of the edge w.r.t. the cell
        switch (edge_orientations[ed_sub_idx]) {
          case lf::mesh::Orientation::positive: {
            for (int j = 0; j < no_int_dof_edge; j++) {
              dofs_[0].push_back(dofs_[1][edge_int_dof_offset + j]);
              cell_dof_offset++;
            }
            break;
          }
          case lf::mesh::Orientation::negative: {
            for (int j = no_int_dof_edge - 1; j >= 0; j--) {
              dofs_[0].push_back(dofs_[1][edge_int_dof_offset + j]);
              cell_dof_offset++;
            }
            break;
          }
        }  // end switch
      }    // end loop over edges

      // Set indices for interior shape functions of the cell
      for (int j = 0; j < no_int_dof_cell; j++) {
        dofs_[0].push_back(dof_idx);
        cell_dof_offset++;
        dof_entities_.push_back(cell);
        dof_idx++;
      }  // end loop over interior dof of cell
    }    // end loop over cells
    // Set sentinel
    offsets_[0][no_cells] = cell_dof_offset;

    // Finally set number of global shape functions
    num_dof_ = dof_idx;
  }  // end constructor

  /**
   * @copydoc DofHandler::GetNoDofs()
   */
  size_type NoDofs() const override { return num_dof_; }

  /**
   * @copydoc DofHandler::GetNoInteriorDofs()
   * @sa InteriorGlobalDofIndices
   */
  size_type NoInteriorDofs(const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GetNoLocalDofs()
   * @sa GlobalDofIndices()
   */
  size_type NoLocalDofs(const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GlobalDofIndices()
   */
  lf::base::RandomAccessRange<const gdof_idx_t> GlobalDofIndices(
      const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::InteriorGlobalDofIndices()
   */
  lf::base::RandomAccessRange<const gdof_idx_t> InteriorGlobalDofIndices(
      const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GetEntity()
   * @sa GlobalDofIndices()
   */
  const lf::mesh::Entity &Entity(gdof_idx_t dofnum) const override {
    LF_VERIFY_MSG(dofnum < dof_entities_.size(),
                  "Illegal dof index " << dofnum << ", max = " << num_dof_);
    return *dof_entities_[dofnum];
  }

  /** @copydoc DofHandler::mesh()
   */
  std::shared_ptr<const lf::mesh::Mesh> Mesh() const override {
    return mesh_p_;
  }

 private:
  /** The mesh on which the degrees of freedom are defined */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  /** The total number of degrees of freedom */
  size_type num_dof_{0};
  /** Vector of entities to which global basis functions are associated */
  std::vector<const lf::mesh::Entity *> dof_entities_;
  /** Internal indexing helper arrays
   *
   * The indices of global shape functions are stored in a long array
   * `dofs_[codim]` for each co-dimension `codim`. The entries in the
   * `offsets_[codim]` array point to the beginning of the dof indices for a
   * particular entity identified via its index in the mesh.
   *
   * @note the length `offsets_[codim]` vectors must exceed the number of
   * corresponding entities by 1.
   */
  /**@{*/
  /** Vector of number of interior dofs for entities */
  std::array<std::vector<size_type>, 3> no_int_dofs_;
  /** Offsets of dof index arrays in `dofs_` vectors */
  std::array<std::vector<size_type>, 3> offsets_;
  /** Vectors of global indices of dofs belonging to entities of different
      co-dimension */
  std::array<std::vector<gdof_idx_t>, 3> dofs_;
  /**@}*/
};

}  // namespace lf::assemble

#endif
