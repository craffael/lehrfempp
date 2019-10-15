/** @file refutils.cc
 * @brief Implementation of debugging utilities
 */

#include "refutils.h"
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <fstream>

namespace lf::refinement {

inline int normalize_idx(glb_idx_t idx) {
  LF_ASSERT_MSG(idx == idx_nil || idx <= std::numeric_limits<int>::max(),
                "Error, trying to convert " << idx << "into an int");
  return (idx == idx_nil) ? int(-1) : int(idx);
}

void WriteMatlabLevel(const MeshHierarchy &hier_mesh, size_type level,
                      std::string filename) {
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
    file << "function [PTPAR,EDPAR,CELLPAR] = " << filename << "()"
         << std::endl;
    file << "% Parent data for a hybid 2D mesh" << std::endl;

    // Mesh on the selected level
    std::shared_ptr<const mesh::Mesh> mesh = hier_mesh.getMesh(level);
    {
      // Output parent information for nodes
      const std::vector<ParentInfo> &pt_parent_info(
          hier_mesh.ParentInfos(level, 2));
      const size_type no_nodes = mesh->NumEntities(2);
      for (int k = 0; k < no_nodes; k++) {
        file << "PTPAR(" << k + 1 << ",:) = ["
             << normalize_idx(pt_parent_info[k].parent_index) << " , "
             << normalize_idx(pt_parent_info[k].child_number) << "];"
             << std::endl;
      }
    }

    {
      // Output parent information for edge
      const std::vector<ParentInfo> &ed_parent_info(
          hier_mesh.ParentInfos(level, 1));
      const size_type no_edges = mesh->NumEntities(1);
      for (int k = 0; k < no_edges; k++) {
        file << "EDPAR(" << k + 1 << ",:) = ["
             << normalize_idx(ed_parent_info[k].parent_index) << " , "
             << normalize_idx(ed_parent_info[k].child_number) << "];"
             << std::endl;
      }
    }

    {
      // Output parent information for cells
      const std::vector<ParentInfo> &cell_parent_info(
          hier_mesh.ParentInfos(level, 0));
      const std::vector<glb_idx_t> &ref_eds(hier_mesh.RefinementEdges(level));
      const size_type no_cells = mesh->NumEntities(0);
      for (int k = 0; k < no_cells; k++) {
        file << "CELLPAR(" << k + 1 << ",:) = ["
             << normalize_idx(cell_parent_info[k].parent_index) << " , "
             << normalize_idx(cell_parent_info[k].child_number) << " , "
             << normalize_idx(ref_eds[k]) << " ];" << std::endl;
      }
    }
  }
}

void WriteMatlab(const MeshHierarchy &hier_mesh, const std::string &basename) {
  const size_type n_levels = hier_mesh.NumLevels();
  for (int level = 0; level < n_levels; level++) {
    // prepare filename
    std::stringstream level_asc;
    level_asc << level;
    std::string filebase = basename + "_L" + level_asc.str();

    // Fetch mesh on the current level
    std::shared_ptr<const mesh::Mesh> mesh = hier_mesh.getMesh(level);

    // Output of mesh data
    io::writeMatlab(*mesh, filebase + ".m");
    // Output of parent/refinement edge information
    WriteMatlabLevel(hier_mesh, level, filebase + "_pi.m");
  }
}

}  // namespace lf::refinement
