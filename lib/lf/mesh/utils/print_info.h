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
 * @brief Diagnostic output operator. Prints info about a mesh.
 * @param &mesh The mesh to print info about
 * @param &o The stream to which this function should output
 */
void PrintInfo(const Mesh &mesh, std::ostream &o);


// Print function for Entity object
/**
 * @brief Diagnostic output operator. Prints info about an entity.
 *
 * @param e The entity to print info about
 * @param stream The stream to which this function should output
 *
 * #### Output levels
 * - Entity::output_ctrl_ == 0: Entity type is printed
 * - Entity::output_ctrl_ > 10: The above and information of codimensions
 * - Entity::output_ctrl_ > 50: The above and information about subentities
 * - Entity::output_ctrl_ > 90: The above and coordinates
 */
void PrintInfo(const lf::mesh::Entity& e, std::ostream& stream);


}  // namespace lf::mesh::utils

#endif  // __a0ec4da7c53444cbb215ff2415c2b3c5
