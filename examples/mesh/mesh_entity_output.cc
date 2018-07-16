/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>

#include "lf/mesh/hybrid2d/entity.h"
#include "lf/mesh/hybrid2dp/triangle.h"
#include "lf/mesh/mesh_factory.h"
#include "lf/mesh/entity.h"
#include "lf/base/base.h"
#include "lf/mesh/tp_triag_mesh_builder.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"



int main() {
  std::cout << "Output of information for mesh entity elements" << std::endl;

  // NOTE TO SELF: There are several test folders in lf/mesh:
  // - mesh/test
  // - mesh/hybrid2d/test
  // - mesh/hybrid2dp/test

  // lf::mesh::hybrid2d::MeshFactory factory(2); // Working, but not entity

  return 0L;
}
