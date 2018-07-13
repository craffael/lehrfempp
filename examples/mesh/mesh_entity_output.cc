/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>

#include "lf/mesh/hybrid2d/entity.h"
#include "lf/mesh/hybrid2dp/triangle.h"
#include "lf/mesh/mesh_factory.h"
#include "lf/mesh/entity.h"
#include "lf/base/base.h"

int main() {
  std::cout << "Output of information for mesh entity elements" << std::endl;

  //auto entity = lf::mesh::hybrid2d::Entity(); // Not working
  auto test = lf::mesh::hybrid2dp::Triangle::RefEl();
  //std::cout << "Print of entity: " << entity << std::endl;

  return 0L;
}
