#include "entity.h"

namespace lf::mesh::hybrid2d {
  // Implementation of sub-entity access methods

  base::RandomAccessRange<const mesh::Entity>
  Edge::SubEntities(char rel_codim) const {
    LF_VERIFY_MSG(rel_codim <= 1,"Edge: rel_codim out of range");
    /*
      Implementation missing
     */
  }

  base::RandomAccessRange<const mesh::Entity>
  Trilateral::SubEntities(char rel_codim) const {
    LF_VERIFY_MSG(rel_codim <= 2,"Trilateral: rel_codim out of range");
    /*
      Implementation missing
     */
  }
  
  base::RandomAccessRange<const mesh::Entity>
  Quadrilateral::SubEntities(char rel_codim) const {
     LF_VERIFY_MSG(rel_codim <= 2,"Quadrilateral: rel_codim out of range");
     /*
      Implementation missing
     */ 
  }

}
