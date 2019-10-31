#ifndef PROJECTS_DPG_PRODUCT_DOFHANDLER
#define PROJECTS_DPG_PRODUCT_DOFHANDLER

/**
 * @file
 * @brief A DOF handler for product spaces.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>

#include "dpg.h"

namespace projects::dpg {

/**
 * @brief see lf::assemble
 */
using dof_map_t = lf::assemble::UniformFEDofHandler::dof_map_t;

/**
 * @brief Dofhandler for uniform finite element spaces that
 * discretizes a cartesian/product space.
 *
 * This class for indices of global shape functions is suitable
 * for situations where the FE-space is based on a Cartesian Product
 * function space and where every geometric entity of a particular
 * type has exactly the same number of shape functions.
 *
 * ###Additional terminology:
 *
 * - We use the term "product space" or "cartesian space" for function spaces
 * \f$ U \f$ that have the following structure \f[ U = U_0 \times U_1 \times
 * \dots \times U_{n-1} \f] where each of the \f$ U_i \f$ is a  function space
 * of scalar  valued functions.
 *
 * - A function in such a space  represents several variables \f$ u = (u_0,
 * \dots u_{n-1}) \f$. We usually use the term component to refer to such a
 * variable. Each component \f$ u_i \f$ is associated to  a unique component
 * index \f$ i \in \{ 0,\dots n-1 \} \f$ representing its position
 *
 *
 * ###Rules for ordering global shape functions.
 *
 * The following rules are used to order global shape functions in such a
 * product space:
 *
 * - Dofs belonging to components of smaller component index are arranged before
 * dofs belonging to components of larger component index.
 * - Within each component the ususal rules for ordering global shape functions
 * specified in lf::assemble::DofHandler apply.
 *
 * ###Rules for ordering local shape functions.
 *
 * - Dofs belonging to components of smaller component index are
 * arranged before dofs associated to components of larger component index.
 * - Within each component the usual rules for ordering local shape functions
 * specified in lf::assemble::DofHandler apply.
 *
 *
 * These rules result in a block structure of both Element matrices and system
 * matrices.
 *
 * This construction of the dof handler
 * also implicitly defines a basis on the Cartesian space based on the bases
 * on the component spaces in the following way: If the dofhandler manages dofs
 * on the Product space \f$ U_h \times  Q_h \f$ , where \f$U_h\f$ and \f$Q_h\f$
 * have bases \f$\{u_i\}\f$ and \f$\{q_j\}\f$, then the basis on the Cartesian
 * space is given by: \f$\{(u_i,0)\} \cup \{(0,q_j)\} \f$, as one would expect.
 */
class ProductUniformFEDofHandler : public lf::assemble::DofHandler {
 public:
  /**
   * @brief Construction from a vector of map objects.
   * @param mesh The underlying mesh
   * @param dofmaps A vector of maps describing the local dof layout
   * for each component. The dofmaps describe the number of interior dofs
   * for every type of entity.
   */
  ProductUniformFEDofHandler(std::shared_ptr<const lf::mesh::Mesh> mesh_p,
                             std::vector<dof_map_t> dofmaps);

  /** @copydoc lf::assemble::DofHandler::NoDofs() */
  [[nodiscard]] size_type NoDofs() const override;

  /**
   * @brief tells the total number of dofs handled, that belong to a
   * component of the Product space.
   * @param component index identifying a component of the product
   * space.
   */
  [[nodiscard]] size_type NoDofs(size_type component) const;

  /** @copydoc lf::assemble::DofHandler::NoLocalDofs() */

  [[nodiscard]] size_type NoLocalDofs(
      const lf::mesh::Entity& entity) const override;

  /**
   * @brief tells the number of dofs subordinate to an entity and
   * belonging to a component of the product space.
   * @param entity entity of the underlying mesh whose number of shape
   * functions is queried
   * @param component  index identifying the component of the product
   * space.
   */
  [[nodiscard]] size_type NoLocalDofs(const lf::mesh::Entity& entity,
                                      size_type component) const;

  /** @copydoc lf::assemble::DofHandler::NoInteriorDofs() */
  [[nodiscard]] size_type NoInteriorDofs(
      const lf::mesh::Entity& entity) const override;

  /**
   * @brief tells the number of dofs associated to an entity and
   * belonging to a component of the product space.
   * @param entity entity of the underlying mesh whose number of shape
   * functions is queried
   * @param component index identifying the component of the product
   * space.
   */
  [[nodiscard]] size_type NoInteriorDofs(const lf::mesh::Entity& entity,
                                         size_type component) const;

  /** @copydoc lf::assemble::DofHandler::GlobalDofIndices() */
  [[nodiscard]] nonstd::span<const gdof_idx_t> GlobalDofIndices(
      const lf::mesh::Entity& entity) const override;

  /** @copydoc lf::assemble::DofHandler::InteriorGlobalDofIndices */
  [[nodiscard]] nonstd::span<const gdof_idx_t> InteriorGlobalDofIndices(
      const lf::mesh::Entity& entity) const override;

  /** @copydoc lf::assemble::DofHandler::Entity() */
  [[nodiscard]] const lf::mesh::Entity& Entity(
      gdof_idx_t dofnum) const override;

  /**
   * @brief retrive the unique index of the component to which a global
   * basis function belongs.
   * @param dofnum global index of the basis function
   *
   * This function returns the unique component of the product space, to which
   * a particular global shape function belongs.
   */
  [[nodiscard]] size_type Component(glb_idx_t dofnum) const;

  /** @copydoc lf::assemble::DofHandler::Mesh() */
  [[nodiscard]] std::shared_ptr<const lf::mesh::Mesh> Mesh() const override;

  /**
   * @brief Returns a DofHandler object, that handles degrees of freedoms
   * belonging to a component.
   * @param component  index of the component for which a dofhandler
   * should be obtained
   *
   * Returns a dofhandler object that manages the global shape functions
   * belonging to a given component. This object then numbers the shape
   * functions of this component from 0 to NoDofs(component)-1.
   *
   * @note The resulting object no longer
   * has any information about the product space from which it was created.
   */
  [[nodiscard]] const lf::assemble::UniformFEDofHandler& ComponentDofHandler(
      size_type component) const {
    return *component_dof_handlers_[component];
  }

  /**
   * @brief Returns the first/smallest index of a local shape functions that is
   * associated/subordinate  and that belongs to a  component on a certain
   * entity.
   * @param entity reference to the entity for which this information is to be
   * fetched.
   * @param component index of the component for which this information is to be
   * fetched.
   *
   * The extended rules for the ordering of local shape functions state that
   * local shape functions belonging to certain component are numbered
   * consecutively. This function returns the smallest local index of a shape
   * function that bleongs to the component and is subordinate/associated to the
   * entity.
   *
   * @note The local shape functions belonging  to  this component and that are
   * subordinate/associated to the entity will have local indieces in the range
   * \f$ \{ \f$ LocalStartIndex(entity,component), \f$ \dots \f$,
   * LocalStartIndex(entity,component) + NoLocalDofs(entity,component) \f$ \}
   * \f$
   */
  [[nodiscard]] ldof_idx_t LocalStartIndex(const lf::mesh::Entity& entity,
                                           size_type component) const;

  /**
   * @brief returns the number of components of the product space.
   */
  [[nodiscard]] size_type NoComponents() const { return num_components_; }

  /**
   * @brief Returns the first/smallest index of a global shape function
   * belonging to a certain component.
   * @param component index of the component for which this information is to be
   * fetched.
   *
   * The extended rules for the ordering of global shape functions state that
   * global shape functions associated with certain component are numbered
   * consecutively. This function returns the smallest global index of a shape
   * functions that belongs to a certain component component.
   *
   * @note The global shape functions associated with this component will have
   * global indices in the range \f$  \{ \f$ Offset(component), \f$\dots\f$ ,
   * Offset(component) + NoDofs(component)\f$ \} \f$
   */
  [[nodiscard]] size_type Offset(size_type component) const {
    return offsets_[component];
  }

 private:
  /**
   * @brief initialization of internal index arrays.
   */
  void initIndexArrays();

  /** mesh on which the degrees of freedom are defined */
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;

  /** dof handlers for each component:*/
  std::vector<std::shared_ptr<lf::assemble::UniformFEDofHandler>>
      component_dof_handlers_;

  /** the total number of components*/
  size_type num_components_;

  /** vector of global indices of dofs subordinate to entities
   * of different topological type*/
  // access via  dof_[codim][entity_idx][local_idx].
  std::array<std::vector<std::vector<gdof_idx_t>>, 3> dofs_;

  /** vectors of global indices of dofs belonging to entities
   * of different topological type*/
  // access via  dof_[codim][entity_idx][local_idx].
  std::array<std::vector<std::vector<gdof_idx_t>>, 3> internal_dofs_;

  /** stores the global offsets for the different components: */
  std::vector<glb_idx_t> offsets_;
};

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_PRODUCT_DOFHANDLER
