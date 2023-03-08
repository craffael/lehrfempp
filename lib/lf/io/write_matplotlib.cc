/**
 * @file
 * @brief Implementation of writeMatplotlib
 * @author Anian Ruoss
 * @date   2018-10-08 18:27:17
 * @copyright MIT License
 */

#include "write_matplotlib.h"

#include <fstream>
#include <iostream>

namespace lf::io {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void writeMatplotlib(const lf::mesh::Mesh &mesh, std::string filename) {
  using dim_t = lf::base::RefEl::dim_t;

  // append .csv to filename if necessary
  if (filename.compare(filename.size() - 4, 4, ".csv") != 0) {
    filename += ".csv";
  }

  std::ofstream file(filename);

  if (file.is_open()) {
    const dim_t dim_mesh = mesh.DimMesh();
    LF_VERIFY_MSG(dim_mesh == 2,
                  "write_matplotlib() only available for 2D meshes");

    // loop through all elements of every codimension
    for (int codim = 0; codim <= dim_mesh; ++codim) {
      for (const lf::mesh::Entity *obj : mesh.Entities(codim)) {
        const size_t obj_idx = mesh.Index(*obj);
        const lf::base::RefEl obj_ref_el = obj->RefEl();
        const lf::geometry::Geometry *obj_geo_ptr = obj->Geometry();
        const Eigen::MatrixXd vertices =
            obj_geo_ptr->Global(obj_ref_el.NodeCoords());

        switch (obj_ref_el) {
          case lf::base::RefEl::kPoint(): {
            file << codim << ',' << obj_idx << ',' << vertices(0, 0) << ','
                 << vertices(1, 0) << std::endl;
            break;
          }
          case lf::base::RefEl::kSegment(): {
            file << codim << ',' << obj_idx << ',';
            // to access points of segment use SubEntities(1)
            for (const auto *sub : obj->SubEntities(codim)) {
              file << mesh.Index(*sub) << ',';
            }
            file << std::endl;

            break;
          }
          case lf::base::RefEl::kTria():
          case lf::base::RefEl::kQuad(): {
            file << codim << ',' << obj_idx << ',';
            // to access points of cell use SubEntities(1)
            for (const auto *sub : obj->SubEntities(codim + 1)) {
              file << mesh.Index(*sub) << ',';
            }
            file << std::endl;

            break;
          }
          default: {
            LF_VERIFY_MSG(false, "Error for object " + std::to_string(obj_idx) +
                                     " in codim " + std::to_string(codim) +
                                     " of type " + obj_ref_el.ToString());
            break;
          }
        }
      }
    }
  }
}

}  // namespace lf::io
