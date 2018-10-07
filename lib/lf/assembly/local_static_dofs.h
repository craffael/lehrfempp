/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Description of location of local shape functions
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef _LF_LOCSTATDOFH_H
#define _LF_LOCSTATDOFH_H

#include <lf/base/base.h>
#include "assembly_types.h"

namespace lf::assemble {

/**
 * @brief Defines local distribution of shape functions = degrees of freedom
 *
 *
 * This _interface_ class should be used to set up dof handlers for uniform
 * finite element spaces, for which the same number of basis functions is
 * associated to every mesh entity of a particular type. The prime example are
 * Lagrangian finite element spaces.
 *
 * ### Example: 2D Quadratic Lagrangian finite elements
 *
 * 1 local dof per point, 1 local dof per edge, 0 local dof for triangles,
 * 1 local dof for quads
 *
 * ### Example: 2D Crouzeix-Raviart element
 *
 * 0 local dof per point, 1 local dof per edge, 0 cell-based dof
 *
 * ### Example: 2D second-order edge elements
 *
 * 0 local dof per node, 2 local dof per edge, 2 local dof per cell
 */
class LocalStaticDOFs {
 public:
  LocalStaticDOFs() = default;
  /** Disabled constructors and assignment operators */
  /**@{*/
  /** @name Disabled */
  LocalStaticDOFs(const LocalStaticDOFs &) = delete;
  LocalStaticDOFs(LocalStaticDOFs &&) = delete;
  LocalStaticDOFs &operator=(const LocalStaticDOFs &) = delete;
  LocalStaticDOFs &operator=(LocalStaticDOFs &&) = delete;
  /**@}*/
  virtual ~LocalStaticDOFs() = default;

  /**
   * @brief Obtain number of local shape functions associated with
   * (sub-)entities of a given type
   *
   * @param refel topological type of the entity.
   *
   * @note the method provides the number of _interior_ local shape functions
   *
   */
  virtual size_type NoLocDofs(lf::base::RefEl refel) const = 0;
  /**
   * @brief The total number of local shape functions belonging to an entity
   * type
   *
   * @param refel topological type of the entity.
   *
   * This method returns the number of all local shape functions whose trace
   * does not vanish on entities of the specified type. This number is equal
   * to the sum of the numbers of interior local shape functions for
   * an entity and all its sub-entities.
   */
  virtual size_type TotalNoLocDofs(lf::base::RefEl refel) const = 0;
  /**
   * @brief topological dimension of mesh
   */
  virtual dim_t Dimension() const = 0;
};

/**
 * @brief Description of uniform distribution of local shape functions for 2D
 * mesh
 *
 * @sa LocalStaticDOFs
 */
class LocalStaticDOFs2D : public LocalStaticDOFs {
 public:
  /**
   * @brief set multiplicities for (sub-)entities of any type
   *
   * @param nodof_point no of local shape functions sitting at nodes
   * @param nodof_segment no of _iterior_ local shape funtions on edges
   * @param nodof_tria  no of _interior_ local shape functions of triangles
   * @param nodof_quad
   */
  LocalStaticDOFs2D(size_type nodof_point, size_type nodof_segment,
                    size_type nodof_tria, size_type nodof_quad)
      : no_loc_dofs_per_refel_() {
    no_loc_dofs_per_refel_[0] = nodof_point;
    no_loc_dofs_per_refel_[1] = nodof_segment;
    no_loc_dofs_per_refel_[2] = nodof_tria;
    no_loc_dofs_per_refel_[3] = nodof_quad;
  }
  ~LocalStaticDOFs2D() override = default;
  LocalStaticDOFs2D(const LocalStaticDOFs2D &) = delete;
  LocalStaticDOFs2D(LocalStaticDOFs2D &&) = delete;
  LocalStaticDOFs2D &operator=(const LocalStaticDOFs2D &) = delete;
  LocalStaticDOFs2D &operator=(LocalStaticDOFs2D &&) = delete;

  /**
   * @copydoc LocalStaticDOFs::NoLocDofs()
   *
   */
  size_type NoLocDofs(lf::base::RefEl refel) const override;
  /**
   * @copydoc LocalStaticDOFs::TotalNoLocDofs()
   */
  size_type TotalNoLocDofs(lf::base::RefEl refel) const override;
  /**
   * @copydoc LocalStaticDOFs::Dimension()
   */
  dim_t Dimension() const override { return 2; }

 private:
  /** number of local shape functions associated with each (sub-)entity
      of a particular co-dimension (multiplicity)*/
  std::array<size_type, 4> no_loc_dofs_per_refel_{};
};  // end class definition

/** @brief Local dof numbers for linear Lagrangian finite elements.
 *
 * For this simplest finite element methods every node owns exactly
 * one global/local shape function.
 */
class LocalLinearLagrangianFE2D : public LocalStaticDOFs2D {
 public:
  LocalLinearLagrangianFE2D() : LocalStaticDOFs2D(1, 0, 0, 0) {}
  LocalLinearLagrangianFE2D(const LocalLinearLagrangianFE2D &) = delete;
  LocalLinearLagrangianFE2D(LocalLinearLagrangianFE2D &&) = delete;
  LocalLinearLagrangianFE2D &operator=(const LocalLinearLagrangianFE2D &) =
      delete;
  LocalLinearLagrangianFE2D &operator=(LocalLinearLagrangianFE2D &&) = delete;
  ~LocalLinearLagrangianFE2D() override = default;
};

}  // namespace lf::assemble

#endif
