#include "lf/base/ref_el.h"
#include <cassert>
#include <iostream>

namespace lf::base {
void foo() {
  //! [convert_to_enum]
  RefElType type = RefEl::kTria();
  //! [convert_to_enum]

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-comparison"
  //! [oneToOneRelation]
  RefElType::kPoint == RefEl::kPoint();
  RefElType::kSegment == RefEl::kSegment();
  RefElType::kTria == RefEl::kTria();
  RefElType::kQuad == RefEl::kQuad();
  //! [oneToOneRelation]
#pragma clang diagnostic pop

  //! [enumConversion]
  RefElType point = RefEl::kPoint();
  auto segment = RefEl(RefElType::kSegment);
  //! [enumConversion]

  //! [nodeCoordStatic]
  // If RefEl not known at compile time:
  auto nodeCoordsDynamic = RefEl::kTria().NodeCoords();

  // If RefEl known at compile time:
  std::vector<Eigen::Vector2d> nodeCoordsCompiletime =
      RefEl::NodeCoords<RefEl::kTria()>();
  //! [nodeCoordStatic]

  //![streamOutput]
  std::cout << RefEl::kSegment();  // prints "kSegment"
  //![streamOutput]
  {
    //![refElUsage]
    // Example usage
    auto triangle = RefEl::kTria();
    assert(triangle.Dimension() == 2);
    assert(triangle.NumNodes() == 3);
    assert(triangle.NodeCoords().col(0) == Eigen::Vector2d(0, 0));
    assert(triangle.NumSubEntities(1) ==
           3);  // Triangle has 3 sub-entities with codim=1 (all segments)

    auto point =
        triangle.SubType(2, 0);  // RefEl of sub-entity with codim=2, index=0
    assert(point == RefEl::kPoint());
    //![refElUsage]
  }
}  // foo
}  // namespace lf::base
