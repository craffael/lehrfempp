/** @file write_matlab.cc */

#include "write_matlab.h"

#include <fstream>

namespace lf::io {

void writeMatlab(const lf::mesh::Mesh& mesh, std::string filename) {
  using size_type = std::size_t;         // unsigned integer
  using dim_t = lf::base::RefEl::dim_t;  // type for dimensions
  const Eigen::MatrixXd zero(Eigen::MatrixXd::Zero(0, 1));

  // Preprocess filename: append .m, if not there already
  {
    const size_type fn_len = filename.size();
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
    file << "function [x,y,TRI,QUAD,EDS] = " << filename << "()" << '\n';
    file << "% Data for an unstructure planar hybrid 2D mesh" << '\n';

    // Obtain topological dimension of the mesh
    const dim_t dim_mesh = mesh.DimMesh();
    LF_VERIFY_MSG(dim_mesh == 2,
                  "write_matlab() only avvailable for 2D meshes");

    // Run through the nodes of the mesh and write x coordinates
    const dim_t node_codim(dim_mesh);
    // Number of nodes of the mesh
    const size_type no_of_nodes = mesh.NumEntities(node_codim);
    file << "x = zeros(" << no_of_nodes << ",1);" << '\n';
    file << "y = zeros(" << no_of_nodes << ",1);" << '\n';

    // Write node coordinates to file
    size_type node_cnt = 0;
    for (const mesh::Entity* node : mesh.Entities(node_codim)) {
      const lf::base::glb_idx_t node_index = mesh.Index(*node);
      const geometry::Geometry* geo_ptr = node->Geometry();
      Eigen::MatrixXd node_coord(geo_ptr->Global(zero));
      file << "x(" << node_index + 1 << ") = " << node_coord(0, 0) << "; "
           << '\n';
      file << "y(" << node_index + 1 << ") = " << node_coord(1, 0) << "; "
           << '\n';
      node_cnt++;
    }
    LF_VERIFY_MSG(node_cnt == no_of_nodes, "Node count mismatch");

    // Write edge information to file
    const size_type no_of_edges = mesh.NumEntities(1);
    file << "EDS = zeros(" << no_of_edges << ",2);" << '\n';
    size_type ed_cnt = 0;
    for (const mesh::Entity* edge : mesh.Entities(1)) {
      const lf::base::glb_idx_t edge_index = mesh.Index(*edge);
      const base::RefEl ref_el = edge->RefEl();
      LF_VERIFY_MSG(ref_el == lf::base::RefEl::kSegment(),
                    "Edge must be a segment");
      // Access endpoints = sub-entities of relative co-dimension 1
      const auto sub_ent = edge->SubEntities(1);
      file << "EDS(" << edge_index + 1 << ",:) = ["
           << mesh.Index(*sub_ent[0]) + 1 << ", " << mesh.Index(*sub_ent[1]) + 1
           << "];" << '\n';
      ed_cnt++;
    }

    // Write cell (entities of co-dimension 0) information to file
    file << "TRI = []; QUAD = [];" << '\n';
    size_type cell_cnt = 0;
    size_type triag_cnt = 0;
    size_type quad_cnt = 0;
    for (const mesh::Entity* e : mesh.Entities(0)) {
      const lf::base::glb_idx_t cell_index = mesh.Index(*e);
      // Access vertices =  sub-entities of relative co-dimension 2
      const auto sub_ent = e->SubEntities(2);
      const base::RefEl ref_el = e->RefEl();
      switch (ref_el) {
        case lf::base::RefEl::kTria(): {
          file << "TRI(" << triag_cnt + 1 << ",:) = ["
               << mesh.Index(*sub_ent[0]) + 1 << ", "
               << mesh.Index(*sub_ent[1]) + 1 << ", "
               << mesh.Index(*sub_ent[2]) + 1 << ", " << cell_index << " ];"
               << '\n';
          triag_cnt++;
          break;
        }
        case lf::base::RefEl::kQuad(): {
          file << "QUAD(" << quad_cnt + 1 << ",:) = ["
               << mesh.Index(*sub_ent[0]) + 1 << ", "
               << mesh.Index(*sub_ent[1]) + 1 << ", "
               << mesh.Index(*sub_ent[2]) + 1 << ", "
               << mesh.Index(*sub_ent[3]) + 1 << ", " << cell_index << " ];"
               << '\n';
          quad_cnt++;
          break;
        }
        default: {
          LF_VERIFY_MSG(false,
                        "write_matlab can only handle triangles and quads");
          break;
        }
      }
      cell_cnt++;
    }
  }
}

}  // namespace lf::io
