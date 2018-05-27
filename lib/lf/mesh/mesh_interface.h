
#ifndef __b86b0c9cb0fd48da931a1f24421b8842
#define __b86b0c9cb0fd48da931a1f24421b8842



#include <lf/base/base.h>

#include "entity.h"

namespace lf::mesh
{
  class Mesh {
  public:

    using size_type = unsigned int;


    /**
     * @brief The dimension of the manifold described by the mesh, or 
     *        equivalently the maximum dimension of the reference elements 
     *        present in the mesh.
     */
    virtual char DimMesh() const = 0;

    
    /**
     * @brief The dimension of the Euclidean space in which the mesh is 
     *        embedded.
     */
    virtual char DimWorld() const = 0;


    /**
     * @brief All entities of a given codimension.
     * @param codim The codimension of the entities that should be returned.
     * @return A base::ForwardRange that can be used to traverse the entities.
     * 
     * @sa Entity
     */
    virtual base::ForwardRange<const Entity> Entities(char codim) const = 0;


    /**
     * @brief The number of Entities that have the given codimension.
     * @param codim The codimension of the entities that should be counted.
     * @return That number of entities that have the given codimension.
     */
    virtual size_type Size(char codim) const = 0;

    // Move this method over to the entity?
    virtual size_type Index(const Entity& e) const = 0;

    
    /**
     * @brief virtual destructor
     */
    virtual ~Mesh() = default;
  };
}


#endif // __b86b0c9cb0fd48da931a1f24421b8842




