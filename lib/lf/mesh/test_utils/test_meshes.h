/** @file test_meshes.h
 * @brief generation routines for test meshes
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef LF_MFT_H
#define LF_MFT_H

#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

/**
 * @brief Utilities for testing sanity of mesh data structures and
 *        tests involving meshes
 */
namespace lf::mesh::test_utils {
/**
 * @brief Generates a simple 2D hybrid test mesh
 *
 * @param selector integer parameter for the selection of test meshes
 * @param scale factor for scaling all the coordinates in the meshes
 *
 * The following line of code provides a pointer to the default test mesh:
 * ~~~
 *   #include "lf/mesh/test_utils/test_meshes.h"
 *   ...
 *   auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
 * ~~~
 *
 * - Test mesh selected with selector = 0: domain \f$[0.3]^2\f$*scale
 * @image html testmesh.png
 * - Test mesh generated when selector = 1;
 * @image html testmesh1.png
 * This mesh contains parallelograms, that is, affine quadrilaterals.
 * - Test mesh generated when selector = 3;
 * @image html testmesh3.png
 * This is a purely triangular mesh
 * - selector = 4: test mesh of [0,3]^2 with triangles and quads
 * - selector = 5: test mesh of [0,3]^2 with triangles and parallelograms, all _affine_ 
 */
std::shared_ptr<lf::mesh::Mesh> GenerateHybrid2DTestMesh(int selector = 0,
                                                         double scale = 1.0);

  static const lf::base::size_type GenerateHybrid2DTestMesh_maxsel = 5;

}  // namespace lf::mesh::test_utils

#endif
