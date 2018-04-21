#include "lf/base/ref_el.h"

using namespace lf::base;

void foo() {
  //! [convert_to_enum]
  RefElType type = RefEl::kTria;
  //! [convert_to_enum]

  //! [oneToOneRelation]
  RefElType::kPoint == RefEl::kPoint;
  RefElType::kSegment == RefEl::kSegment;
  RefElType::kTria == RefEl::kTria;
  RefElType::kQuad == RefEl::kQuad;
  //! [oneToOneRelation]

  //! [enumConversion]
  RefElType point = RefEl::kPoint;
  RefEl segment = RefElType::kSegment;
  //! [enumConversion]

  //! [nodeCoordStatic]
  // If RefEl not known at compile time:
  std::vector<Eigen::VectorXd> nodeCoordsDynamic = RefEl::kTria.NodeCoords();

  // If RefEl known at compile time:
  std::vector<Eigen::Vector2d> nodeCoordsCompiletime = RefEl::NodeCoords<RefEl::kTria>();
  //! [nodeCoordStatic]
}
