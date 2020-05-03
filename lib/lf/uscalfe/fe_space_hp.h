#ifndef LF_USCALFE_FE_SPACE_HP_H_
#define LF_USCALFE_FE_SPACE_HP_H_

#include "hp_fe.h"


namespace lf::uscalfe {

/**
 * @brief Lagrangian Finite Element Space of arbitrary degree
 */
template <typename SCALAR>
class FeSpaceHP : public lf::uscalfe::UniformScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  FeSpaceHP() = delete;
  FeSpaceHP(const FeSpaceHP &) = delete;
  FeSpaceHP(FeSpaceHP &&) noexcept = default;
  FeSpaceHP &operator=(const FeSpaceHP &) = delete;
  FeSpaceHP &operator=(FeSpaceHP &&) noexcept = default;

  /**
   * @brief Constructor: Sets up the dof handler
   * @param mesh_p A shared pointer to the underlying mesh (immutable)
   * @param N The polynomial degree of the Finite Element Space
   */
  explicit FeSpaceHP(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, unsigned N)
      : lf::uscalfe::UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<FeHPTria<SCALAR>>(N),
            std::make_shared<FeHPQuad<SCALAR>>(N),
            std::make_shared<FeHPSegment<SCALAR>>(N),
            std::make_shared<lf::uscalfe::FeLagrangePoint<SCALAR>>(N)) {
	// This constructor does very unorthodox things, but there is no easy way around it
	// with the current state of LehrFEM++
	// Permute the edge dofs such that they have a global instead of a local ordering
	const auto& dofh = lf::uscalfe::UniformScalarFESpace<SCALAR>::LocGlobMap();
	for (auto cell : mesh_p->Entities(0)) {
	    const auto dofidxs = dofh.GlobalDofIndices(*cell);
	    // Extreamly evil but necessary as UniformFEDofHandler does not expose its dofs_ array
	    lf::assemble::gdof_idx_t* edge_dofidx = const_cast<lf::assemble::gdof_idx_t*>(dofidxs.data());
	    edge_dofidx += cell->RefEl().NumSubEntities(2);
	    const auto orient = cell->RelativeOrientations();
	    const auto edges = cell->SubEntities(1);
	    // Iterate over all edges and reverse their dofs if orient[i] is negative
	    for (int i = 0 ; i < edges.size() ; ++i) {
		const auto num_edge_dofs = dofh.NumInteriorDofs(*(edges[i]));
		if (orient[i] == lf::mesh::Orientation::negative) {
		    std::reverse(edge_dofidx, edge_dofidx+num_edge_dofs);
		}
		edge_dofidx += num_edge_dofs;
	    }
	}
    }

  ~FeSpaceHP() override = default;
};

}   // end namespace lf::uscalfe


#endif // LF_USCALFE_FE_SPACE_HP_H_
