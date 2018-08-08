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
void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename);


// enum
enum TikzFlags {
    EdgeNumbering = 1,
    CellNumbering = 2,
    NodeNumbering = 4
};

/*
// https://stackoverflow.com/questions/1448396/how-to-use-enums-as-flags-in-c
inline TikzFlags operator|(TikzFlags a, TikzFlags b){
    return static_cast<TikzFlags>(static_cast<int>(a) | static_cast<int>(b));
}
*/



} // namespace lf::mesh::utils



#endif // WRITE_TIKZ_H


