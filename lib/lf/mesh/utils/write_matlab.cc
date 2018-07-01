/** @file write_matlab.cc */

#include "write_matlab.h"

#include <fstream>

namespace lf::mesh::utils {

void writeMatlab(const lf::mesh::Mesh &mesh, std::string filename) {
  using size_type = std::size_t;         // unsigned integer
  using dim_t = lf::base::RefEl::dim_t;  // type for dimensions
  const Eigen::MatrixXd zero(Eigen::MatrixXd::Zero(0, 1));

  // Preprocess filename: append .m, if not there already
  {
    size_type fn_len = filename.size();
    if ((fn_len > 1) && (filename[fn_len - 2] != '.') &&
        (filename[fn_len - 1] != 'm')) {
      filename += ".m";
    }
  }
  // Open file for writing
  std::ofstream file(filename);
  if (file.is_open()) {
    // Start regular writing operations
    // Remove trailing .m
    filename.pop_back();
    filename.pop_back();
    file << "function [x,y,TRI.EDS] = " << filename << "()" << std::endl;

    // Obtain topological dimension of the mesh
    const dim_t dim_mesh = mesh.DimMesh();
    LF_VERIFY_MSG(dim_mesh == 2,
                  "write_matlab() only avvailable for 2D meshes");

    // Run through the nodes of the mesh and write x coordinates
    const dim_t node_codim(dim_mesh);
    // Number of nodes of the mesh
    const size_type no_of_nodes = mesh.Size(node_codim);
    file << "x = zeros(" << no_of_nodes << ",1)" << std::endl;
    file << "y = zeros(" << no_of_nodes << ",1)" << std::endl;

    // Write node coordinates to file
    size_type node_cnt = 0;
    for (const Entity &node : mesh.Entities(node_codim)) {
      const geometry::Geometry *geo_ptr = node.Geometry();
      Eigen::MatrixXd node_coord(geo_ptr->Global(zero));
      file << "x(" << node_cnt << ") = " << node_coord(0, 0) << "; "
           << std::endl;
      file << "y(" << node_cnt << ") = " << node_coord(1, 0) << "; "
           << std::endl;
      node_cnt++;
    }
    LF_VERIFY_MSG(node_cnt == no_of_nodes, "Node count mismatch");

    // Write edge information to file
    const size_type no_of_edges = mesh.Size(1);
    file << "EDS = zeros(" << no_of_edges << ",2);";
    size_type ed_cnt = 0;
    for (const Entity &edge : mesh.Entities(1)) {
      base::RefEl ref_el = edge.RefEl();
      LF_VERIFY_MSG(ref_el == lf::base::RefEl::kSegment(),
                    "Edge must be a segment");
      // Access endpoints = sub-entities of relative co-dimension 1
      const auto sub_ent = edge.SubEntities(1);
      file << "EDS(" << ed_cnt << ",:) = [" << mesh.Index(sub_ent[0]) << ", "
           << mesh.Index(sub_ent[1]) << "];" << std::endl;
      ed_cnt++;
    }

    // Write cell (entities of co-dimension 0) information to file
    const size_type no_of_cells = mesh.Size(0);
    file << "TRI = zeros(" << no_of_cells << ",3);" << std::endl;
    size_type cell_cnt = 0;
    for (const Entity &triangle : mesh.Entities(0)) {
      base::RefEl ref_el = triangle.RefEl();
      LF_VERIFY_MSG(ref_el == lf::base::RefEl::kTria(),
                    "wite_matlab can handle triangular meshes only!");
      // Access vertices =  sub-entities of relative co-dimension 2
      const auto sub_ent = triangle.SubEntities(2);
      file << "TRI(" << cell_cnt << ",:) = [" << mesh.Index(sub_ent[0]) << ", "
           << mesh.Index(sub_ent[1]) << ", " << mesh.Index(sub_ent[2]) << "];"
           << std::endl;
      cell_cnt++;
    }
  }
}

}  // namespace lf::mesh::utils
