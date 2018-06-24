#ifndef HYBRID2D_P_H_
#define HYBRID2D_P_H_

/**
 * @brief An alternative implementation of a hybrid2d mesh manager that uses
 *        Pointers to store sub-entity relations.
 *
 * At the moment this mesh manager doesn't work yet completly but it is destined
 * to replace the `lf::mesh::hybrid2d` mesh manager.
 */
namespace lf::mesh::hybrid2dp {}

#include "mesh.h"
#include "mesh_factory.h"
#include "point.h"
#include "quad.h"
#include "segment.h"
#include "triangle.h"

#endif  // HYBRID2D_P_H_
