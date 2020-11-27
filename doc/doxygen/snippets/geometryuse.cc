/**
 * @file
 * @brief Doxygen snippets to show usage of LehrFEM++'s geometry functionality
 * @author Ralf Hiptmair
 * @date Jan 28, 2020
 * @copyright MIT License
 */

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <iostream>

// The following snippet is copied from lecturedemomesh.cc
//! [corners]
void PrintGeometryInfo(const lf::mesh::Mesh &mesh, lf::base::dim_t codim) {
  // loop over all entities of the specified codimension
  for (const lf::mesh::Entity *ent : mesh.Entities(codim)) {
    // Number of nodes = number of corner points
    const lf::base::size_type num_nodes = ent->RefEl().NumNodes();
    // Obtain pointer to geometry object associated with entity
    const lf::geometry::Geometry *geo_ptr = ent->Geometry();
    // Fetch coordinates of corner points in packed format \cref{par:coords}
    Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
    std::cout << ent->RefEl() << "(" << mesh.Index(*ent) << ") pts: ";
    for (int l = 0; l < num_nodes; ++l) {
      std::cout << l << " =[" << corners.col(l).transpose() << "], ";
    }
    std::cout << std::endl;
  }
}
//! [corners]
