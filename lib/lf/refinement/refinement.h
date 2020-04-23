#ifndef _LF_REFINEMENT_H_
#define _LF_REFINEMENT_H_

/**
 * @file refinement.h
 * @brief Definitions of data structures related to refinement of the mesh
 *
 */

#include <lf/mesh/mesh.h>
#include "hybrid2d_refinement_pattern.h"
#include "mesh_hierarchy.h"
#include "refutils.h"

/**
 * @brief tools for regular or local refinement of 2D hybrid meshes
 *
 * LehrFEM++ supports the _local_ refinement of hybrid meshes ensuring mesh
 * conformity in the process. Edge marking is used to single out cells that must
 * be refined. Other cells may also undergo refinement in order to avoid hanging
 * nodes.
 */
namespace lf::refinement {}  // namespace lf::refinement

#endif
