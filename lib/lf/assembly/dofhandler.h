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
   * @brief A general (interface) class for DOF handling
   *
   */
  class DofHandler {
  public:
    /** 
     * @brief The constructor just stores the mesh pointer
     */
    DofHandler(std::shared_ptr<lf::mesh::Mesh> mesh):mesh_(std::move(mesh)) {
      LF_VERIFY_MSG(mesh_ != nullptr,"Invalid mesh pointer");
    }
    virtual ~DofHandler(void) = default;
  protected:
    DofHandler(const DofHandler &) = default;
    DofHandler(DofHandler &&) = default;
    DofHandler &operator=(const DofHandler &) = default;
    DofHandler &operator=(DofHandler &&) = default;
  public:
    /** 
     * @brief total number of dof's handled by the object
     */
    size_type GetNoDofs(void) const { return num_dof_; }
    
    /** 
     * @brief access to indices of global dof's belonging to an entity 
     * 
     * @param entity reference to the entity for which the dof's are to be 
     *        fetched. This entity must belong to the underlying mesh.
     *
     * The basis functions of every finite element space must be associated with a 
     * unique geometric entity. Conversely, every entity can possess a finite number
     * of basis functions = degrees of freedom. This functions return the global indices
     * of all basis functions associated with the entity _and its sub-entitites_.
     */
    virtual std::vector<gdof_idx_t>
    GetGlobalDofs(const lf::mesh::Entity &entity) const = 0;

    /**
     * @brief tells the number of degrees of freedom subordinated to an entity
     * 
     * @param entity reference to an entity of the underlying mesh
     *
     * This method informs about the lenght of the vector returned by GetGlobalDofs(). 
     * @sa GetGlobalDofs()
     */
    virtual size_type GetNoDofs(const lf::mesh::Entity &entity) const = 0;
    
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
    virtual const lf::mesh::Entity &GetEntity(gdof_idx_t dofnum) const =  0;
    
  protected:
    /** The mesh on which the degrees of freedom are defined */
    std::shared_ptr<lf::mesh::Mesh> mesh_;
    /** The total number of degrees of freedom */
    size_type num_dof_{0};
  };
  

} // end namespace

#endif
