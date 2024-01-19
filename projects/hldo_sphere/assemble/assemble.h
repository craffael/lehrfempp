#ifndef HLDO_SPHERE_ASSEMBLE_H
#define HLDO_SPHERE_ASSEMBLE_H

#include <laplace_matrix_provider.h>
#include <load_vector_provider.h>
#include <mass_matrix_provider.h>
#include <rot_whitney_one_div_matrix_provider.h>
#include <whitney_one_curl_curl_matrix_provider.h>
#include <whitney_one_grad_matrix_provider.h>
#include <whitney_one_mass_matrix_provider.h>
#include <whitney_one_vector_provider.h>
#include <whitney_two_mass_matrix_provider.h>
#include <whitney_two_vector_provider.h>

namespace projects::hldo_sphere {
/**
 * @brief Collection of matrix and vector element providers
 *
 * Containing all necessary matrix and vector providers for the
 * Hodge Laplacians and the Dirac operator
 *
 * Details for the mathematical derivations can be found in the thesis (linked
 * in the README.md) chapter 4.
 */

namespace assemble {}
}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_ASSEMBLE_H
