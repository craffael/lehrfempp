/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of local computations for FE for
 * 2nd-order linear elliptic boundary value problems.
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "loc_comp_ellbvp.h"

namespace lf::fe {

std::shared_ptr<spdlog::logger> diffusion_element_matrix_provider_logger =
    base::InitLogger("lf::fe::reaction_element_matrix_provider_logger");

std::shared_ptr<spdlog::logger> mass_element_matrix_provider_logger =
    base::InitLogger("lf::fe::mass_element_matrix_provider_logger");

std::shared_ptr<spdlog::logger> mass_edge_matrix_provider_logger =
    base::InitLogger("lf::fe::mass_edge_matrix_provider_logger");

std::shared_ptr<spdlog::logger> scalar_load_element_vector_provider_logger =
    base::InitLogger("lf::fe::scalar_load_element_vector_provider_logger");

std::shared_ptr<spdlog::logger> scalar_load_edge_vector_provider_logger =
    base::InitLogger("lf::fe::scalar_load_edge_vector_provider_logger");

}  // end namespace lf::fe
