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
 * @brief Diagnostic output operator. Prints info about a mesh.
 * @param &mesh The mesh to print info about
 * @param &o The stream to which this function should output
 *
 * #### Output levels
 * - > 10: also output entity information
 * - > 90: also output subentity information
 */
void PrintInfo(std::ostream& o, const lf::mesh::Mesh& mesh, int ctrl = 11);

// Print function for Entity object
/**
 * @brief Prints info about an entity.
 *
 * @param e The entity to print info about
 * @param stream The stream to which this function should output
 *
 * #### Output levels
 * - output_ctrl == 0: Entity type is printed
 * - output_ctrl > 10: The above and information of codimensions
 * - output_ctrl > 50: The above and information about subentities
 * - output_ctrl > 90: The above and coordinates
 */
void PrintInfo(std::ostream& stream, const lf::mesh::Entity& e,
               int output_ctrl = 0);

}  // namespace lf::mesh::utils

#endif  // __a0ec4da7c53444cbb215ff2415c2b3c5
