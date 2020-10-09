/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief All implementation in header file
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include "loc_comp_norms.h"

namespace lf::uscalfe {

std::shared_ptr<spdlog::logger> mesh_function_l2_norm_difference_logger =
    base::InitLogger("lf::uscalfe::mesh_function_l2_norm_difference_logger");

std::shared_ptr<spdlog::logger> mesh_function_l2_gradient_difference_logger =
    base::InitLogger(
        "lf::uscalfe::mesh_function_l2_gradient_difference_logger");

}  // namespace lf::uscalfe
