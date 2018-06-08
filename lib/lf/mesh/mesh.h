/** @file mesh.h
 * @brief Header file for the mesh topology facilities of LehrFEM++
 */
#ifndef __7fd0772dc0d1474492ad46dc2b8cebbb
#define __7fd0772dc0d1474492ad46dc2b8cebbb

#include "entity.h"
#include "mesh_builder.h"
#include "mesh_interface.h"

namespace lf::mesh
{
  /** 
   * @brief Function for testing mesh indexing
   * @param mesh A reference to the mesh to be checked
   * @return False, if there is a problem with indexing
   * This function tests whether all entities of a mesh are indexed 
   * consecutively. 
   */
  bool testEntityIndexing(const Mesh &mesh);
}

#endif  // __7fd0772dc0d1474492ad46dc2b8cebbb
