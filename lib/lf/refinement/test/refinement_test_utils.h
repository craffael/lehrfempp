/**
 * @file
 * @brief Check that the father child relations that are returned from the
 *        refinement module are correct.
 * @author Raffael Casagrande
 * @date   2018-08-12 01:00:35
 * @copyright MIT License
 */

#ifndef LF_REFTEST_UTILS_H
#define LF_REFTEST_UTILS_H

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

namespace lf::refinement::test {
void checkFatherChildRelations(const MeshHierarchy &mh,
                               base::size_type father_level);

void checkGeometryInParent(const MeshHierarchy &mh,
                           base::size_type father_level);

}  // namespace lf::refinement::test

#endif
