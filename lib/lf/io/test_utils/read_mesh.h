/**
 * @file
 * @brief Defines utility functions to read a test mesh from a file
 *        that is in the git repository
 * @author Raffael Casagrande
 * @date   2018-07-21 11:56:52
 * @copyright MIT License
 */

#ifndef __1bd576875bb04c168d82c67cc451cbe7
#define __1bd576875bb04c168d82c67cc451cbe7
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <string>
#include "lf/io/gmsh_reader.h"

namespace lf::io::test_utils {

/// Retrieves the name of a mesh that lies in the lib/lf/io/test/msh_files/
/// folder
std::string getMeshPath(std::string mesh_name);

GmshReader getGmshReader(std::string mesh_name, base::dim_t dim_world);

}  // namespace lf::io::test_utils

#endif  // __1bd576875bb04c168d82c67cc451cbe7
