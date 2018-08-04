/** @file refutils.h
 * @brief utility functions for debugging refinement methods
 */

#include "refinement_hierarchy.h"

#ifndef LF_REFUTILS_H_H
#define LF_REFUTILS_H_H

namespace lf::refinement {

/**
 * @brief Generate MATLAB function providing parent/child information
 * @param hier_mesh reference to a valid MeshHierarchy object
 * @param level refinement level to be output
 * @filename output file in which to store the MATLAB code
 */
void WriteMatlabLevel(const MeshHierarchy &hier_mesh, size_type level,
                      std::string filename);

/**
 * @brief Generate MATLAB code describing a multilevel hierarchy
 *        of meshes
 * @param hier_mesh reference to a valid MeshHierarchy object
 * @filename basename for output file, of which there will be _several_
 */
void WriteMatlab(const MeshHierarchy &hier_mesh,const std::string &basename);

}  // namespace lf::refinement

#endif
