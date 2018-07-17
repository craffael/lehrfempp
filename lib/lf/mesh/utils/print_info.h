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
 * @brief diagnostic output operator
 */
void PrintInfo(const Mesh &mesh, std::ostream &o);

/*
// Print function for Entity object
void PrintInfoEntity(const lf::mesh::Entity&, std::ostream& stream);
// ?? void PrintInfoEntity(const lf::mesh::Entity& e, std::ostream &stream);
*/



}  // namespace lf::mesh::utils

#endif  // __a0ec4da7c53444cbb215ff2415c2b3c5
