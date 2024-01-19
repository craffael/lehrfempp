#ifndef INCGea25e7ba269040a1876408a8c67712b5
#define INCGea25e7ba269040a1876408a8c67712b5

#include <lf/mesh/mesh.h>

namespace lf::mesh::test_utils {

/**
 * @brief Makes sure that the coordinates of the nodes obtained through
 *        `e.SubEntities()[i].Geometry()` match the ones obtained through
 *         `e.Geometry()`
 * @param e The entity to check.
 */
void checkGeometryOrientation(const Entity& e);

}  // namespace lf::mesh::test_utils

#endif  // INCGea25e7ba269040a1876408a8c67712b5
