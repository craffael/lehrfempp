#ifndef LF_LD_REF_H
#define LF_LD_REF_H
/**
 * @file
 * @brief Simple LehrFEM demo code for refinement
 * @author Ralf Hiptmair
 * @date   March 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/refinement/refinement.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lecturedemo {

/** @brief Creates a hierarchy of meshes by regular refinement
 *
 * @param mesh reference to a mesh object
 *
 * This function just creates a mesh hierarchy by regular refinement
 * and outputs the numbers of various entities on different levels.
 * Does not serve any other purpose.
 */
void regrefMeshSequence(const std::shared_ptr<lf::mesh::Mesh>& mesh_p,
                        int refsteps);

/** @brief Demonstration of regular refinement: driver function
 */
void lecturedemorefine();

}  // namespace lecturedemo

#endif
