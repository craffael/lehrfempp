#ifndef __37e385afbd3b4b1dba8611fb71787822
#define __37e385afbd3b4b1dba8611fb71787822
#include "../base/forward_range.h"
#include "geometry.h"
#include "../base/ref_el.h"


namespace lf::mesh {

class Entity {

public:

  virtual char Codim() = 0;

  /**
   * @brief Return all sub entities of this entity that have the given 
   *        codimension (w.r.t. this entity!)
   * @param codim The codim w.r.t. this entity
   * @return 
   */
  virtual base::ForwardRange<Entity> SubEntities(char codim) = 0;

  /**
   * @brief Describes the geometry of this entity.
   * @return A point to a Geometry object that will remain valid for as long
   *         as the Mesh remains valid.
   */
  virtual Geometry* Geometry() = 0;

  /**
   * @brief Describes the reference element type of this entity.
   * @return An object of type base::RefEl.
   */
  virtual base::RefEl RefEl() = 0;

};


}

#endif // __37e385afbd3b4b1dba8611fb71787822
