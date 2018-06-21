#ifndef __2633e7bf8dd548de839f74075687e81x
#define __2633e7bf8dd548de839f74075687e81x

#include <lf/mesh/mesh.h>

namespace lf::mesh::test_utils {
/**
 * @brief Function for testing consistency of subentities
 * @param e a (reference) to a mesh entity
 *
 * This function verifies the consistency of the local topology, that is
 * whether subentities of sub-entities are also subentities of the entity.
 * 
 */
void checkLocalTopology(const Entity& e);

}  // namespace lf::mesh::test_utils

#endif  // __2633e7bf8dd548de839f74075687e81x
