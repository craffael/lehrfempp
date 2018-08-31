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

#ifndef _LF_DOFHD_H
#define _LF_DOFHD_H

#include "local_static_dofs.h"

namespace lf::assemble {
/*
 * @brief A general (_interface_) class for DOF handling
 *
 */
class DofHandler {
 public:
  /**
   * @brief The constructor just stores the mesh pointer
   */
  DofHandler(void) {}
  virtual ~DofHandler(void) = default;

  DofHandler(const DofHandler &) = delete;
  DofHandler(DofHandler &&) = delete;
  DofHandler &operator=(const DofHandler &) = delete;
  DofHandler &operator=(DofHandler &&) = delete;
  /**
   * @brief total number of dof's handled by the object
   */
  virtual size_type GetNoDofs(void) const = 0;

  /**
   * @brief access to indices of global dof's belonging to an entity
   *
   * @param ref_el_type type of entity for which to retrieve dof indices
   * @param entity_index unique mesh index for an entity of the specified type
   *
   * The basis functions of every finite element space must be associated with a
   * unique geometric entity. Conversely, every entity can possess a finite
   * number of basis functions = degrees of freedom. This functions return the
   * global indices of all basis functions associated with the entity _and its
   * sub-entitites_.
   */
  virtual lf::base::RandomAccessRange<const gdof_idx_t> GetGlobalDofs(
      lf::base::RefEl ref_el_type, glb_idx_t entity_index) const = 0;

  /**
   * @brief access to indices of global dof's belonging to an entity
   *
   * @param entity reference to the entity for which the dof's are to be
   *        fetched. This entity must belong to the underlying mesh.
   *
   * @sa GetGlobalDofs()
   */
  virtual lf::base::RandomAccessRange<const gdof_idx_t> GetGlobalDofs(
      const lf::mesh::Entity &entity) const = 0;

  /**
   * @brief tells the number of degrees of freedom subordinated to an entity
   *
   * @param entity reference to an entity of the underlying mesh
   *
   * This method informs about the lenght of the vector returned by
   * GetGlobalDofs().
   * @sa GetGlobalDofs()
   */
  virtual size_type GetNoLocalDofs(const lf::mesh::Entity &entity) const = 0;

  /** @brief returns number of shape functions covering an entity
   *  @sa getNoLocalDofs()
   */
  virtual size_type GetNoLocalDofs(lf::base::RefEl ref_el_type,
                                   glb_idx_t entity_index) const = 0;

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
   * @sa GetGlobalDofs()
   */
  virtual const lf::mesh::Entity &GetEntity(gdof_idx_t dofnum) const = 0;

  /** @brief Acess to underlying mesh object
   *
   * Every DofHandler object has to be associated with a unique mesh.
   * All entities passed to the DofHandler object must belong to that mesh.
   */
  virtual const lf::mesh::Mesh &getMesh(void) const = 0;
};

/**
 * @brief Dofhandler for uniform finite element spaces
 *
 * This management class for indices of global shape functions
 * is suitable for situations where every geometric entity of a particular
 * type has exactly the same number of shape functions belonging to it.
 */
class UniformFEDofHandler : public DofHandler {
 public:
  /** @brief Initialization of global index arrays
   */
  UniformFEDofHandler(std::shared_ptr<lf::mesh::Mesh> mesh,
                      const LocalStaticDOFs &locdof);

  /**
   * @copydoc DofHandler::GetNoDofs()
   */
  virtual size_type GetNoDofs(void) const override { return num_dof_; }

  /**
   * @copydoc DofHandler::GetGlobalDofs()
   */
  virtual lf::base::RandomAccessRange<const gdof_idx_t> GetGlobalDofs(
      lf::base::RefEl ref_el_type, glb_idx_t entity_index) const override;

  virtual lf::base::RandomAccessRange<const gdof_idx_t> GetGlobalDofs(
      const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GetNoLocalDofs()
   * @sa GetGlobalDofs()
   */
  virtual size_type GetNoLocalDofs(
      const lf::mesh::Entity &entity) const override;

  /**
   * @copydoc DofHandler::GetNoLocalDofs()
   * @sa GetGlobalDofs()
   */
  virtual size_type GetNoLocalDofs(lf::base::RefEl ref_el_type,
                                   glb_idx_t entity_index) const override;

  /**
   * @copydoc DofHandler::GetEntity()
   * @sa GetGlobalDofs()
   */
  virtual const lf::mesh::Entity &GetEntity(gdof_idx_t dofnum) const override {
    LF_VERIFY_MSG(dofnum < dof_entities_.size(),
                  "Illegal dof index " << dofnum << ", max = " << num_dof_);
    return *dof_entities_[dofnum];
  }

  /** @copydoc DofHandler::getMesh()
   */
  virtual const lf::mesh::Mesh &getMesh(void) const override { return *mesh_; }

 private:
  const size_type kNodeOrd = 0;
  const size_type kEdgeOrd = 1;
  const size_type kCellOrd = 2;

  /** The mesh on which the degrees of freedom are defined */
  std::shared_ptr<lf::mesh::Mesh> mesh_;
  /** The total number of degrees of freedom */
  size_type num_dof_{0};
  /** Vector of entities to which global basis functions are associated */
  std::vector<const lf::mesh::Entity *> dof_entities_;
  /** Vectors of global indices of dofs belonging to entities of different
      topological type */
  std::array<std::vector<gdof_idx_t>, 3> dofs_;
  /** Number of dofs belonging to entities of a particular type */
  std::array<size_type, 3> no_dofs_;
  /** Number of shape functions covering different cell types */
  size_type num_dofs_tria_{0}, num_dofs_quad_{0};
};

}  // namespace lf::assemble

#endif
