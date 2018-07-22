#ifndef WRITE_TIKZ_H
#define WRITE_TIKZ_H

#include <lf/mesh/mesh.h>
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"


#include <string>

namespace lf::mesh::utils {


// Write doxygen documentation here
/**
 * @brief Writes entity data to file in LaTeX/TikZ format.
 *
 * @param entity the entity to be stored to file
 * @param filename name of output file
 *
 * This function writes a file of code, which included in LaTeX draws your entity using TikZ Graphics.
 *
 *
 */

//Declare
//void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename);

void writeTikZ(const lf::mesh::Entity &e, std::string filename);




} // namespace lf::mesh::utils



#endif // WRITE_TIKZ_H


