/** @file mesh/mesh.h
 * @brief Header file for the mesh topology facilities of LehrFEM++
 */
#ifndef __7fd0772dc0d1474492ad46dc2b8cebbb
#define __7fd0772dc0d1474492ad46dc2b8cebbb

/**
 * @brief Defines a set of interface classes that define a mesh manager and
 *        provides mesh-related tools that build on these interfaces.
 *
 * @see lf::mesh::hybrid2d for a concrete implementation of a %mesh manager.
 */
namespace lf::mesh {}

#include "entity.h"
#include "mesh_factory.h"
#include "mesh_interface.h"
#include "tp_triag_mesh_builder.h"

#endif  // __7fd0772dc0d1474492ad46dc2b8cebbb
