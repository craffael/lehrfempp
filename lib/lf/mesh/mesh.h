/** @file mesh/mesh.h
 * @brief Header file for the mesh topology facilities of LehrFEM++
 */
#ifndef INCG7fd0772dc0d1474492ad46dc2b8cebbb
#define INCG7fd0772dc0d1474492ad46dc2b8cebbb

/**
 * @brief Defines a set of interface classes that define a mesh manager and
 *        provides mesh-related tools that build on these interfaces.
 *
 * @see lf::mesh::hybrid2d for a concrete implementation of a %mesh manager.
 *
 * For a discussion of the functionality and design of mesh data structures see
 * [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{sec:meshdata}. A mathematical definition of a finite element mesh is
 * given in [Lecture
 * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
 * @lref{sec:meshes}.
 */
namespace lf::mesh {}

#include "entity.h"
#include "mesh_factory.h"
#include "mesh_interface.h"

#endif  // INCG7fd0772dc0d1474492ad46dc2b8cebbb
