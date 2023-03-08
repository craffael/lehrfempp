/**
 * @file
 * @brief Defines utility functions to read a test mesh from a file
 *        that is in the git repository
 * @author Raffael Casagrande
 * @date   2018-07-21 11:56:52
 * @copyright MIT License
 */

#ifndef INCG1bd576875bb04c168d82c67cc451cbe7
#define INCG1bd576875bb04c168d82c67cc451cbe7
#include <lf/mesh/hybrid2d/hybrid2d.h>

#include <string>

#include "lf/io/gmsh_reader.h"

namespace lf::io::test_utils {

/**
 * @brief Retrieve the full path to the file
 * `lib/lf/io/test/msh_files/<mesh_name>`
 */
std::string getMeshPath(std::string mesh_name);

/**
 * @brief Get a GmshReader from the file `lib/lf/io/test/msh_files/<mesh_name>`.
 * @param mesh_name The name of the mesh file.
 * @param dim_world The world dimension of the mesh.
 */
GmshReader getGmshReader(std::string mesh_name, base::dim_t dim_world);

}  // namespace lf::io::test_utils

#endif  // INCG1bd576875bb04c168d82c67cc451cbe7
