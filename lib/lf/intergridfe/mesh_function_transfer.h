#ifndef LF_INTERGRIDFE_MESH_FUNCTION_TRANSFER_H
#define LF_INTERGRIDFE_MESH_FUNCTION_TRANSFER_H

#include <Eigen/Dense>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/uscalfe.h>



namespace lf::intergridfe {

    template<typename SCALAR_FE_COARSE, typename SCALAR_FE_FINE, typename SCALAR_COEFF>
    [[nodiscard]] Eigen::Matrix<SCALAR_FE_FINE, Eigen::Dynamic, 1>
    transferCoarseFineFE(const lf::refinement::MeshHierarchy &mh,
	    const lf::uscalfe::UniformScalarFESpace<SCALAR_FE_FINE> & fes_fine,
	    const lf::uscalfe::MeshFunctionFE<SCALAR_FE_COARSE, SCALAR_COEFF> &mf) {
	// Find the index of the coarse mesh in mh
	const lf::base::size_type mesh_index = [&]() {
	    lf::base::size_type idx = 0;
	    for (const auto mesh : mh.getMeshes()) {
		if (mesh == mf.getMesh()) {
		    return idx;
		}
		++idx;
	    }
	    return idx;
	}();
	// Assert whether the user provided a valid mesh
	LF_ASSERT_MSG(mesh_index < mh.NumLevels()-1, "Invalid Mesh provided to function.");
	// Get the refinement info of the meshes
	const auto& child_info = mh.CellChildInfos(mesh_index);
	// Iterate over all coarse cells and lift the coefficients
	const auto mesh_fine = fes_fine.Mesh();
	const auto dofh_fine = fes_fine.LocGlobMap();
	const lf::uscalfe::FeSpaceLagrangeO1<SCALAR_FE_FINE> fes_lagr_o1_fine(mesh_fine);   // Used for coordinate transformations
	Eigen::Matrix<SCALAR_FE_FINE, Eigen::Dynamic, 1> coeffs_fine(dofh_fine.NumDofs());
	const auto cells = mf.getMesh().Entities(0);
	for (size_t cell_idx = 0 ; cell_idx < cells.size() ; ++cell_idx) {
	    const auto& cell = *cells[cell_idx];
	    // Reconstruct the refinement pattern
	    const lf::refinement::Hybrid2DRefinementPattern refpat(cell.RefEl(), child_info[cell_idx].ref_pat_, child_info[cell_idx].anchor_);
	    // Iterate over the child cells and set their dof coefficients
	    const auto cp_cells = refpat.ChildPolygons(0);
	    for (lf::base::size_type child_idx = 0 ; child_idx < cp_cells.size() ; ++child_idx) {
		const auto& child_cell = *(mesh_fine->EntityByIndex(child_info[cell_idx].child_cell_idx[child_idx]));
		const auto scalar_ref_fel = fes_fine.ShapeFunctionLayout(child_cell.RefEl());
		const auto sfl_o1 = fes_lagr_o1_fine.ShapeFunctionLayout(child_cell);
		const auto eval_nodes = (cp_cells[child_idx] * sfl_o1.EvalReferenceShapeFunctions(sfl_o1.EvaluationNodes()) / refpat.LatticeConst()).eval();
		const auto eval_node_values = mf(cell, eval_nodes);
		// Convert the vector output of the mesh function to a eigen matrix again
		using mf_t = lf::uscalfe::MeshFunctionFE<SCALAR_FE_COARSE, SCALAR_COEFF>;
		using mf_ret_t = lf::uscalfe::MeshFunctionReturnType<mf_t>;
		const Eigen::Map<Eigen::Matrix<mf_ret_t, 1, Eigen::Dynamic>> nodal_values(eval_node_values.data(), eval_node_values.size());
		const auto dofs = scalar_ref_fel.NodalValuesToDofs(nodal_values.template cast<SCALAR_FE_FINE>());
		const auto glob_dof_idxs = dofh_fine.GlobalDofIndices(child_cell);
		for (lf::base::size_type dof_idx = 0 ; dof_idx < dofs.size() ; ++dof_idx) {
		    coeffs_fine[glob_dof_idxs[dof_idx]] = dofs[dof_idx];
		}
	    }
	}
	// Return the dofvector on the fine mesh
	return coeffs_fine;
    }

}   // end namespace intergridfe



#endif // LF_INTERGRIDFE_MESH_FUNCTION_TRANSFER_H
