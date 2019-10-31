#ifndef PROJECTS_DPG_DPG_H
#define PROJECTS_DPG_DPG_H
#include <lf/uscalfe/uscalfe.h>

/**
 * @brief Contains functionality for the implementation of DPG methods
 */
namespace projects::dpg {

/** @brief  Type for indices into global matrices/vectors */
using gdof_idx_t = lf::uscalfe::gdof_idx_t;

/** @brief  Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::uscalfe::ldof_idx_t;

/** @brief Type for vector length/matrix sizes */
using size_type = lf::uscalfe::size_type;

/** @brief Tpe for (co)-dimensions */
using dim_t = lf::uscalfe::dim_t;

/** @brief Type for global index of entities */
using glb_idx_t = lf::uscalfe::glb_idx_t;

/** @brief Type for indexing of sub-entities*/
using sub_idx_t = lf::uscalfe::sub_idx_t;

}  // namespace projects::dpg

// shape functions
//#include "discontinuous_fe_constant.h"
//#include "lagr_fe_quadratic.h"
//#include "lagr_fe_cubic.h"

// shape function decorators
//#include "discontinuous_scalar_reference_finite_element.h"
//#include "trace_scalar_reference_finite_element.h"

// DPG local quanties
//#include "dpg_element_matrix_provider.h"
//#include "dpg_element_vector_provider.h"

// Sub element matrix providers
//#include "loc_comp_dpg.h"

// product spaces
//#include "product_dofhandler.h"

#endif  // PROJECTS_DPG_DPG_H
