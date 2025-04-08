#ifndef HYBRID2D_P_H_
#define HYBRID2D_P_H_

/**
 * @brief A namespace for managing 2D hybrid meshes using pointers to store
 *        sub-entity relations.
 *
 * This namespace provides an implementation for handling 2D hybrid meshes,
 * which include various types of elements such as points, segments, triangles,
 * and quadrilaterals.
 */
namespace lf::mesh::hybrid2d {}

#include "mesh.h"
#include "mesh_factory.h"
#include "point.h"
#include "quad.h"
#include "segment.h"
#include "triangle.h"

#endif  // HYBRID2D_P_H_
