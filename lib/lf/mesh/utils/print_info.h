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
 * - Entity::output_ctrl_ == 0: Derived type of entity and type of entity is printed
 * - Entity::output_ctrl_ > 0: The above and geometry of the entity is printed
 * - Entity::output_ctrl_ > 10: The above and geometry of subentities is printed
 */
void PrintInfo(const lf::mesh::Entity& e, std::ostream& stream);

/*
// void PrintInfo(const lf::mesh::Entity& e, const lf::mesh::Mesh& m, std::ostream &stream);
*/



}  // namespace lf::mesh::utils

#endif  // __a0ec4da7c53444cbb215ff2415c2b3c5
