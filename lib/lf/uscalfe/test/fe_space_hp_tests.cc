/**
 * @file
 * @brief Check that the FeSpaceHP class works as expected
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>


namespace lf::uscalfe::test {

TEST(fe_space_hp, continuity) {
    for (int selector = 0 ; selector <= 8 ; ++selector) {
	// Get a hybrid test mesh on [0, 3]^2
	const auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
	for (unsigned p = 1 ; p <= 20 ; ++p) {
	    // Generate a FeSpaceHP with order p on the mesh
	    const auto fe_space = std::make_unique<lf::uscalfe::FeSpaceHP<double>>(mesh, p);
	    // Test the continuity of each basis function associated with an edge
	    for (auto&& cell : mesh->Entities(0)) {
		const auto sfl_cell = fe_space->ShapeFunctionLayout(*cell);
		const auto cell_refel = cell->RefEl();
		const auto cell_num_nodes = cell_refel.NumNodes();
		const auto cell_nodes = cell_refel.NodeCoords();
		const auto edges = cell->SubEntities(1);
		const auto orient = cell->RelativeOrientations();
		for (int i = 0 ; i < edges.size() ; ++i) {
		    const Eigen::RowVectorXd edge_eval_coords = Eigen::RowVectorXd::LinSpaced(p+1, 0, 1);
		    const Eigen::MatrixXd cell_eval_coords = [&]() {
			Eigen::MatrixXd result(2, p+1);
			const Eigen::Vector2d node0 = cell_nodes.col((i+0)%cell_num_nodes);
			const Eigen::Vector2d node1 = cell_nodes.col((i+1)%cell_num_nodes);
			result.row(0) = Eigen::RowVectorXd::Constant(p+1, node0[0]) + (node1[0]-node0[0])*edge_eval_coords;
			result.row(1) = Eigen::RowVectorXd::Constant(p+1, node0[1]) + (node1[1]-node0[1])*edge_eval_coords;
			return result;
		    }();
		    const auto edge = edges[i];
		    const auto sfl_edge = fe_space->ShapeFunctionLayout(*edge);
		    const auto rsf_edge = sfl_edge->EvalReferenceShapeFunctions(edge_eval_coords);
		    const auto rsf_cell = sfl_cell->EvalReferenceShapeFunctions(cell_eval_coords);
		    // Compare each basis function on the edge
		    ASSERT_TRUE(sfl_edge->NumRefShapeFunctions(0) == sfl_cell->NumRefShapeFunctions(1)) << "selector=" << selector << " p=" << p << " cell=" << mesh->Index(*cell) << " edge=" << i << std::endl;
		    for (int rsf_idx = 0 ; rsf_idx < sfl_edge->NumRefShapeFunctions(0) ; ++rsf_idx) {
			const Eigen::RowVectorXd rsf_edge_eval = rsf_edge.row(2+rsf_idx);
			Eigen::RowVectorXd rsf_cell_eval;
			if (orient[i] == lf::mesh::Orientation::positive) {
			    rsf_cell_eval = rsf_cell.row(cell_num_nodes+(p-1)*i+rsf_idx);
			}
			else {
			    rsf_cell_eval = rsf_cell.row(cell_num_nodes+(p-1)*(i+1)-rsf_idx-1).reverse();
			}
			const double max_diff = (rsf_edge_eval - rsf_cell_eval).array().abs().maxCoeff();
			ASSERT_TRUE(max_diff < 1e-10) << "selector=" << selector << " p=" << p << " cell=" << mesh->Index(*cell) << " edge=" << i << " rsf_idx=" << rsf_idx << "\nrsf_edge_eval=[" << rsf_edge_eval << "]\nrsf_cell_eval=[" << rsf_cell_eval << "]" << std::endl;
		    }
		}
	    }
	}
    }
}

}   // end namespace lf::uscalfe::test
