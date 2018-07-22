/** @file refutils.cc
 * @brief Implementation of debugging utilities
 */

#include "refutils.h"

#include <fstream>

namespace lf::refinement {

  inline glb_idx_t normalize_idx(glb_idx_t idx) {
    return (idx == idx_nil)?-1:idx;
  }
  
  void WriteMatlabLevel(const MeshHierarchy &hier_mesh,
			size_type level,std::string filename) {
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
      file << "function [PTPAR,EDPAR,CELLPAR] = " << filename << "()" << std::endl;
      file << "% Parent data for a hybid 2D mesh" << std::endl;

      // Mesh on the selected level
      const lf::mesh::Mesh &mesh(hier_mesh.getMesh(level));

      {
	// Output parent information for nodes
	const std::vector<ParentInfo> &pt_parent_info(hier_mesh.parent_info(level,2));
	const size_type no_nodes = mesh.Size(2);
	for (int k=0; k < no_nodes; k++) {
	  file << "PTPAR(" << k+1 << ",:) = ["
	       << normalize_idx(pt_parent_info[k].parent_index_) << " , "
	       << normalize_idx(pt_parent_info[k].child_number_) << "];" << std::endl;
	}
      }

      {
	// Output parent information for edge
	const std::vector<ParentInfo> &ed_parent_info(hier_mesh.parent_info(level,1));
	const size_type no_edges = mesh.Size(1);
	for (int k=0; k < no_edges; k++) {
	  file << "EDPAR(" << k+1 << ",:) = ["
	       << normalize_idx(ed_parent_info[k].parent_index_) << " , "
	       << normalize_idx(ed_parent_info[k].child_number_) << "];" << std::endl;
	}
      }

      {
	// Output parent information for cells
	const std::vector<ParentInfo> &cell_parent_info(hier_mesh.parent_info(level,0));
	const std::vector<glb_idx_t> &ref_eds(hier_mesh.refinement_edges(level));
	const size_type no_cells = mesh.Size(0);
	for (int k=0; k < no_cells; k++) {
	  file << "CELLPAR(" << k+1 << ",:) = ["
	       << normalize_idx(cell_parent_info[k].parent_index_) << " , "
	       << normalize_idx(cell_parent_info[k].child_number_)
	       << normalize_idx(ref_eds[k]) << " ];" << std::endl;
	}
      }
    }
  }
    

  void WriteMatlab(const MeshHierarchy &hier_mesh,
		   std::string filename) {

    }
  
  } // end namespace
