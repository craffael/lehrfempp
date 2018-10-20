/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions for the initialization of special data sets
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef _LF_SPEC_DATA_
#define _LF_SPEC_DATA_

#include "codim_mesh_data_set.h"

namespace lf::mesh::utils {

/**
 * @brief store number of adjacent super-entities
 *
 * @param mesh_p reference to underlying mesh
 * @param codim_sub co-dimension of the queried entities
 * @param codim_super _relative_ co-dimension (with positive sign) of
 *                      super entities.
 *
 * For each entity of a given co-dimension, this function counts the number of
 * adjacent super-entities of some smaller co-dimension.
 *
 * ### Example
 *
 * If, for a 2D mesh we want to count the number of cells adjacent to edges,
 * we have to specify codim_sub = 1, codim_super = 1!
 *
 * @note codim_super is _relative_ with flipped sign!
 */

CodimMeshDataSet<lf::base::size_type> countNoSuperEntities(
    const std::shared_ptr<const Mesh>& mesh_p, lf::base::dim_t codim_sub,
    lf::base::dim_t codim_super);

}  // namespace lf::mesh::utils

#endif
