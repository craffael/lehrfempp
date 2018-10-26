/** @file test_meshes.h
 * @brief generation routines for test meshes
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef LF_MFT_H
#define LF_MFT_H

#include <lf/mesh/hybrid2dp/hybrid2dp.h>
#include <lf/mesh/mesh.h>

/**
 * @brief Utilities for testing sanity of mesh data structures and
 *        tests involving meshes
 */
namespace lf::mesh::test_utils {
/**
 * @brief Generates a simple 2D hybrid test mesh
 *
 * @param selector integer parameter for the selection of test meshes
 *
 * The following line of code provides a pointer to the test mesh:
 * ~~~
 *   #include "lf/mesh/test_utils/test_meshes.h"
 *   ...
 *   auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
 * ~~~
 *
 * - Test mesh selected with selector = 0:
 * @image html testmesh.png
 * - Test mesh generated when selector = 1;
 * @image html testmesh1.png
 * This mesh contains parallelograms, that is, affine quadrilaterals.
 * - Test mesh generated when selector = 3;
 * @image html testmesh3.png
 * This is a purely triangular mesh
 */
std::shared_ptr<lf::mesh::Mesh> GenerateHybrid2DTestMesh(int selector = 0);

}  // namespace lf::mesh::test_utils

#endif
