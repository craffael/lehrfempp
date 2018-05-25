
#ifndef __b86b0c9cb0fd48da931a1f24421b8842
#define __b86b0c9cb0fd48da931a1f24421b8842



#include <lf/base/base.h>

#include "entity.h"

namespace lf::mesh
{
  class Mesh {
  public:

    using size_t = int;

    virtual char DimMesh() const = 0;

    virtual char DimWorld() const = 0;

    virtual base::ForwardRange<Entity> Entities(char codim) const = 0;

    virtual size_t Size(char codim) const = 0;

    // Move this method over to the entity?
    virtual size_t Index(const Entity& e) const = 0;

    virtual ~Mesh() {}
  };
}


#endif // __b86b0c9cb0fd48da931a1f24421b8842




