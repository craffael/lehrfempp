#include "entity.h"
#include "lf/mesh/utils/utils.h"

namespace lf::mesh {
// Introduce output_ctrl_

CONTROLDECLARECOMMENT(Entity, output_ctrl_, "output_ctrl_",
                      "Diagnostics control for Mesh/Entity");

std::ostream& operator<<(std::ostream& stream, const lf::mesh::Entity& entity) {
  if (Entity::output_ctrl_ == 0) {
    return stream << entity.RefEl();
  } else {
    lf::mesh::utils::PrintInfo(entity, stream);
    return stream;
  }
}  // end output  operator <<

}  // namespace lf::mesh
