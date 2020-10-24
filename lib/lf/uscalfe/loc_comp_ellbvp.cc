/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of local computations for Lagrange FE for
 * 2nd-order linear elliptic boundary value problems.
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "loc_comp_ellbvp.h"

namespace lf::uscalfe {

std::shared_ptr<spdlog::logger>
    reaction_diffusion_element_matrix_provider_logger = base::InitLogger(
        "lf::uscalfe::reaction_diffusion_element_matrix_provider_logger");

std::shared_ptr<spdlog::logger> mass_edge_matrix_provider_logger =
    base::InitLogger("lf::uscalfe::mass_edge_matrix_provider_logger");

std::shared_ptr<spdlog::logger> scalar_load_element_vector_provider_logger =
    base::InitLogger("lf::uscalfe::scalar_load_element_vector_provider_logger");

std::shared_ptr<spdlog::logger> scalar_load_edge_vector_provider_logger =
    base::InitLogger("scalar_load_edge_vector_provider_logger");

}  // end namespace lf::uscalfe
