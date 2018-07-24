/**
 * @file
 * @brief Defines a number of PrintInfo functions that output mesh objects
 *        to streams
 * @author Raffael Casagrande
 * @date   2018-07-01 01:32:05
 * @copyright MIT License
 */

#ifndef __a0ec4da7c53444cbb215ff2415c2b3c5
#define __a0ec4da7c53444cbb215ff2415c2b3c5

#include <lf/mesh/mesh.h>

namespace lf::mesh::utils {

  /**
   * @brief output control variable for mesh output
   *
   * - > 10: also output entity information
   * - > 90: also output subentity information
   *
   * @note: this static variable is intialized to the value 100
   */   
  extern int printinfo_ctrl;
  
/**
 * @brief diagnostic output operator
 */
void PrintInfo(const Mesh &mesh, std::ostream &o);
}  // namespace lf::mesh::utils

#endif  // __a0ec4da7c53444cbb215ff2415c2b3c5
