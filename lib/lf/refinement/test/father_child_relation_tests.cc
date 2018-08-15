/**
 * @file
 * @brief Check that the father child relations that are returned from the
 *        refinement module are correct.
 * @author Raffael Casagrande
 * @date   2018-08-12 01:00:35
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/test_utils/read_mesh.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>

namespace lf::refinement::test {

void checkFatherChildRelations(const MeshHierarchy& mh,
                               base::size_type father_level) {
  // Obtain pointers to two consecutive meshes
  std::shared_ptr<const mesh::Mesh> father_mesh = mh.getMesh(father_level);
  std::shared_ptr<const mesh::Mesh> child_mesh = mh.getMesh(father_level + 1);

  // Array with information about offsprings of a point
  const std::vector<lf::refinement::PointChildInfo> &
    point_child_infos = mh.PointChildInfos(father_level);
  // Array with information about children of edges
  const std::vector<lf::refinement::EdgeChildInfo> &
    edge_child_infos = mh.EdgeChildInfos(father_level);
  // Array with information about child entities for cells
  const std::vector<lf::refinement::CellChildInfo> &
    cell_child_infos = mh.CellChildInfos(father_level);

  // Array with information about parents of nodes (co-dimension = 2)
  const std::vector<lf::refinement::ParentInfo> &
    point_father_infos = mh.ParentInfos(father_level + 1, 2);

  // ----------------------------------------
  // I: check parent-child relations for nodes
  // ----------------------------------------

  auto origin = Eigen::VectorXd::Zero(0);
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<bool>>
    child_found =  lf::mesh::utils::make_CodimMeshDataSet<bool>(father_mesh, 2, false);
  for (const lf::mesh::Entity &cp : child_mesh->Entities(2)) {
    const lf::base::glb_idx_t node_index = child_mesh->Index(cp);
    // Info record on parent of current node
    const lf::refinement::ParentInfo&
      father_info = point_father_infos[node_index];
    // Pointer to potential parent; every node of a child mesh must have a parent
    const lf::mesh::Entity *fp = father_info.parent_ptr;
    EXPECT_NE(fp,nullptr)
      << "Mising parent for node " << node_index << " of child mesh";
    // Index of father
    const lf::base::glb_idx_t father_index = father_mesh->Index(*fp);
    EXPECT_EQ(father_index, father_info.parent_index)
      << "Index-pointer mismatch for node " << node_index;
    
    switch (fp->Codim()) {
    case 2: { // father entity is a point
      EXPECT_EQ(father_info.child_number, 0)
	<< "Node can only have a single child node";
      EXPECT_TRUE(cp.Geometry()->Global(origin).isApprox(
          fp->Geometry()->Global(origin)))
	<< "Father and child node at different locations!";
      // Check whether current node is registered with its parent
      const lf::refinement::PointChildInfo &pci(point_child_infos[father_index]);
      EXPECT_EQ(pci.child_point_idx,node_index)
	<< "Father/child node index mismatch: "
	<< pci.child_point_idx << " <-> " << node_index;
      break;
    }
    case 1: {
      // father entity is a segment ->calculate distance to segment:
      // Note: works for straight edges only
      auto a = fp->Geometry()->Global(Eigen::VectorXd::Zero(1));
      auto b = fp->Geometry()->Global(Eigen::VectorXd::Ones(1));
      auto c = cp.Geometry()->Global(origin);
      double dist = std::abs((Eigen::MatrixXd(2, 2) << (b - a), (c - a))
			     .finished().determinant())/(b - a).norm();
      EXPECT_LT(dist, 1e-10)
	<< "Child point not located on straight segment";

      // Check whether current node is registered with its parent edge
      const lf::refinement::EdgeChildInfo &eci(edge_child_infos[father_index]);
      const std::vector<lf::base::glb_idx_t> &cpi(eci.child_point_idx);
      EXPECT_NE(std::find(cpi.begin(),cpi.end(),node_index),cpi.end())
	<< "Child node " << node_index << " not registered with parent edge";
      break;
    }
    case 0: {
      // Parent entity is a cell
      // TODO: Test whether point is located inside the cell
      // (for planar cells with straight edges)

      // Check whether current node is registered with its parent cell
      const lf::refinement::CellChildInfo &cci(cell_child_infos[father_index]);
      const std::vector<lf::base::glb_idx_t> &cpi(cci.child_point_idx);
      EXPECT_NE(std::find(cpi.begin(),cpi.end(),node_index),cpi.end())
	<< "Child node " << node_index << " not registered with parent cell";
      break;
    }
    default: {
      FAIL() << "Illegal co-dimension " << fp->Codim();
      break;
    }
    } // end switch
  } // end loop over nodes of fine mesh

  // ----------------------------------------
  // II: check parent-child relations for edges and cells
  // ----------------------------------------
  // // Array with information about parents of edges (co-dimension = 1)
  // const std::vector<lf::refinement::ParentInfo> &
  //   edge_father_infos = mh.ParentInfos(father_level + 1, 1);
  // // Array with information about parents of cells (co-dimension = 0)
  // const std::vector<lf::refinement::ParentInfo> &
  //   cell_father_infos = mh.ParentInfos(father_level + 1, 0);

  // Run through all edges/cells of the fine mesh
  for (lf::base::dim_t codim = 0; codim <= 1; ++codim) {
    for (const lf::mesh::Entity &ent : child_mesh->Entities(codim)) {
      const lf::base::glb_idx_t entity_index = child_mesh->Index(ent);
      // Obtain info record about parent of current entity
      const lf::refinement::ParentInfo &entity_parent_info =
	(mh.ParentInfos(father_level+1,codim)).at(entity_index);
      // One field gives a pointer to the parent entity
      const lf::mesh::Entity *entity_parent_ptr = entity_parent_info.parent_ptr;
      EXPECT_NE(entity_parent_ptr,nullptr)
	<< "Invalid pointer to parent of entity " << entity_index;
      if (entity_parent_ptr != nullptr) {
	// The other field tells the index of the parent in the parent mesh
	const lf::base::glb_idx_t entity_parent_index = entity_parent_info.parent_index;
	// Check consistency of index of father entity
	EXPECT_EQ(entity_parent_index,father_mesh->Index(*entity_parent_ptr))
	  << "Inconsistent parent index: " << entity_parent_index << " <-> "
	  << father_mesh->Index(*entity_parent_ptr);
	const lf::base::RefEl parent_type = entity_parent_ptr->RefEl();
	EXPECT_GE(codim,entity_parent_ptr->Codim())
	  << "Parent co-dimenension " << (int)entity_parent_ptr->Codim()
	  << " < " << (int)codim;
	// The parent can be either an edge or a cell
	switch (parent_type) {
	case lf::base::RefEl::kPoint(): {
	  ADD_FAILURE() << "A point cannot be the parent of an edge";
	  break;
	}
	case lf::base::RefEl::kSegment(): {
	  // Parent is another edge
	  const lf::refinement::EdgeChildInfo &
	    eci(edge_child_infos.at(entity_parent_index));
	  const std::vector<lf::base::glb_idx_t> &cei(eci.child_edge_idx);
	  EXPECT_NE(std::find(cei.begin(),cei.end(),entity_index),cei.end())
	    << "Child entity " << entity_index << " not registered with edge";
	  break;
	}
	case lf::base::RefEl::kTria():
	case lf::base::RefEl::kQuad(): {
	  // Parent is a cell
	  const lf::refinement::CellChildInfo &
	    cci(cell_child_infos.at(entity_parent_index));
	  switch (codim) {
	  case 1: {
	    // Child is an edge
	    const std::vector<lf::base::glb_idx_t> &cei(cci.child_edge_idx);
	    EXPECT_NE(std::find(cei.begin(),cei.end(),entity_index),cei.end())
	      << "Child edge " << entity_index << " not registered with parent cell";
	    break;
	  }
	  case 0: {
	    // Child is a cell as well the parent
	    const std::vector<lf::base::glb_idx_t> &ccidx(cci.child_cell_idx);
	    EXPECT_NE(std::find(ccidx.begin(),ccidx.end(),entity_index),ccidx.end())
	      << "Child cell " << entity_index << " not registered with parent cell";
	    break;
	  }
	  } // end switch codim
	  break;
	}
	default: {
	  FAIL() << "Unknown entity type";
	  break;
	}
	} // end switch
      }
    } // end loop over entities
  } // end loop over codimensions
}

TEST(lf_refinement, FatherChildRelations) {
  auto gmsh_reader =
      io::test_utils::getGmshReader("two_element_hybrid_2d.msh", 2);
  auto base_mesh = gmsh_reader.mesh();

  MeshHierarchy mh(base_mesh,
                   std::make_shared<mesh::hybrid2dp::MeshFactory>(2));

  // mark all edges of the triangle for refinement:
  auto marks = mesh::utils::make_CodimMeshDataSet(base_mesh, 1, false);
  auto triangle = std::find_if(
      base_mesh->Entities(0).begin(), base_mesh->Entities(0).end(),
      [](const auto& e) { return e.RefEl() == base::RefEl::kTria(); });

  for (auto& edge : triangle->SubEntities(1)) {
    (*marks)(edge) = true;
  }

  mh.MarkEdges([&](auto& mesh, auto& e) { return (*marks)(e); });
  mh.RefineMarked();

  {
    // For debug purposes: print the meshs into a vtk file
    io::VtkWriter writer0(base_mesh, "mesh0.vtk");
    io::VtkWriter writer1(mh.getMesh(1), "mesh1.vtk");
  }

  auto child_mesh1 = mh.getMesh(1);
  checkFatherChildRelations(mh, 0);
}

}  // namespace lf::refinement::test
