/**
 * @file
 * @brief Driver function for simple LehrFEM++ demo
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/utils/lambda_mesh_data_set.h"

namespace lecturedemo {
using size_type = lf::base::size_type;
using dim_t = lf::base::dim_t;
using glb_idx_t = lf::base::glb_idx_t;
using sub_idx_t = lf::base::sub_idx_t;

/** @brief traverses entities of specified co-dimension and
 *         prints information about them.
 * @param mesh reference to a read-only LehrFEM++ mesh
 * @param codim co-dimension of interest.
 */
  int traverseEntities(const lf::mesh::Mesh &mesh, dim_t codim);
  
/** @brief function printing the topological relationships of mesh entities
 *  @param mesh read-only LehrFEM++ mesh object
 * @param codim co-dimension of interest.
 */
  void scanTopology(const lf::mesh::Mesh &mesh, dim_t codim);
  
};  // namespace lecturedemo
