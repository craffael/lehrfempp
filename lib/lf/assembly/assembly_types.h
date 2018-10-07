/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Some typedefs used by the LehrFEM++ assembly facilities
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef _LF_ASSEMTYPES_H
#define _LF_ASSEMTYPES_H

#include <lf/base/base.h>

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

}  // namespace lf::assemble

#endif
