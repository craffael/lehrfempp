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
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>

namespace lf::assemble {
/** Type for indices into global matrices/vectors */
using gdof_idx_t = Eigen::Index;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = Eigen::Index;
/** Type for vector length/matrix sizes */
using size_type = lf::base::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::base::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::base::glb_idx_t;

/**
 * @brief Defines local distribution of shape functions = degrees of freedom
 *
 *
 * This class should be used to set up dof handlers for uniform finite element
 * spaces, for which the same number of basis functions is associated to every
 * mesh entity of a particular type. The prime example are Lagrangian finite
 * element spaces.
 *
 * ###Example: 2D Quadratic Lagrangian finite elements
 *
 * 1 local dof per point, 1 local dof per edge, 0 local dof for triangles,
 * 1 local dof for quads
 */
class LocalStaticDOFs {
 public:
  virtual ~LocalStaticDOFs(void) = default;
  LocalStaticDOFs(void);
  LocalStaticDOFs(const LocalStaticDOFs &) = delete;
  LocalStaticDOFs(LocalStaticDOFs &&) = delete;
  LocalStaticDOFs &operator=(const LocalStaticDOFs &) = delete;
  LocalStaticDOFs &operator=(LocalStaticDOFs &&) = delete;

  /**
   * @brief Obtain number of local shape functions associated with
   * (sub-)entities of a given type
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
   * does not vanish on entities of the specified type.
   */
  virtual size_type TotalNoLocDofs(lf::base::RefEl refel) const = 0;
  /**
   * @brief dimension of mesh
   */
  virtual dim_t Dimension(void) const = 0;
};

/**
 * @brief Description of uniform distribution of local shape functions for 2D
 * mesh
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
                    size_type nodof_tria, size_type nodof_quad) {
    no_loc_dofs_per_refel_[0] = nodof_point;
    no_loc_dofs_per_refel_[1] = nodof_segment;
    no_loc_dofs_per_refel_[2] = nodof_tria;
    no_loc_dofs_per_refel_[3] = nodof_quad;
  }
  virtual ~LocalStaticDOFs2D(void) = default;
  LocalStaticDOFs2D(const LocalStaticDOFs2D &) = delete;
  LocalStaticDOFs2D(LocalStaticDOFs2D &&) = delete;
  LocalStaticDOFs2D &operator=(const LocalStaticDOFs2D &) = delete;
  LocalStaticDOFs2D &operator=(LocalStaticDOFs2D &&) = delete;

  /**
   * @copydoc LocalStaticDOFs::NoLocDofs()
   *
   */
  virtual size_type NoLocDofs(lf::base::RefEl refel) const override;
  /**
   * @copydoc LocalStaticDOFs::TotalNoLocDofs()
   */
  virtual size_type TotalNoLocDofs(lf::base::RefEl refel) const override;
  /**
   * @copydoc LocalStaticDOFs::Dimension()
   */
  virtual dim_t Dimension(void) const override { return 2; }

 private:
  /** number of local shape functions associated with each (sub-)entity
      of a particular co-dimension (multiplicity)*/
  std::array<size_type, 4> no_loc_dofs_per_refel_;
};  // end class definition

class LocalLinearLagrangianFE2D : public LocalStaticDOFs2D {
 public:
  LocalLinearLagrangianFE2D(void) : LocalStaticDOFs2D(1, 0, 0, 0) {}
};

}  // namespace lf::assemble

#endif
