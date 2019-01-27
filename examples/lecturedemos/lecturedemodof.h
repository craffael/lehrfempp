#ifndef LF_LD_DOF_H
#define LF_LD_DOF_H
/**
 * @file
 * @brief Functions for simple LehrFEM++ demos + sample codes
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include "lecturedemo.h"
#include "lf/assemble/assemble.h"
#include "lf/mesh/utils/utils.h"
#include "lf/io/io.h"

namespace lecturedemo {

/**
 * @brief output of information stored in DofHandler
 * @parm dofh reference to a DofHandler object
 * @sa std::ostream &lf::assmeble::operator<<(std::ostream &o, const DofHandler
 * &dof_handler)
 */
void printDofInfo(const lf::assemble::DofHandler &dofh);

/**
 * @brief Driver routine for demos for LehrFEM++ assemble module
 */
void lecturedemodof();

}  // namespace lecturedemo

#endif  // LF_LD_DOF_H
