/** @file test_meshes.h
 * @brief generation routines for test meshes
 */

#ifndef LF_MFT_H
#define LF_MFT_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2dp/hybrid2dp.h>

namespace lf::mesh::test_utils {
  /**
   * @brief Generates a simple 2D hybrid test mesh
   */
  
  std::shared_ptr<lf::mesh::Mesh> GenerateHybrid2DTestMesh(void);

}

#endif
