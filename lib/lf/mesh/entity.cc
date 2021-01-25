#include "entity.h"

#include "lf/mesh/utils/utils.h"

namespace lf::mesh {

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
  return stream << entity.RefEl();
}  // end output  operator <<

}  // namespace lf::mesh
