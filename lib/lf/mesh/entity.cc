#include "entity.h"
#include "lf/mesh/utils/utils.h"

namespace lf::mesh {
// Introduce output_ctrl_

CONTROLDECLARECOMMENT(Entity, output_ctrl_, "output_ctrl_",
                      "Diagnostics control for Mesh/Entity");

int to_sign(Orientation o) { return static_cast<int>(o); }
char to_char(Orientation o) {
  switch (o) {
    case lf::mesh::Orientation::positive: {
      return '+';
    }
    case lf::mesh::Orientation::negative: {
      return '-';
    }
  }
  return 0;
}

std::ostream& operator<<(std::ostream& stream, const lf::mesh::Entity& entity) {
  if (Entity::output_ctrl_ == 0) {
    return stream << entity.RefEl();
  }
  utils::PrintInfo(entity, stream);
  return stream;
}  // end output  operator <<

}  // namespace lf::mesh
