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
 * @brief diagnostic output operator. Prints info about a mesh.
 */
void PrintInfo(const Mesh &mesh, std::ostream &o);



// Print function for Entity object
// lf::mesh::Entity?
void PrintInfo(const Entity& e, std::ostream& stream);
/*
// void PrintInfo(const lf::mesh::Entity& e, const lf::mesh::Mesh& m, std::ostream &stream);
*/



}  // namespace lf::mesh::utils

#endif  // __a0ec4da7c53444cbb215ff2415c2b3c5
