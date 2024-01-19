#ifndef HLDO_SPHERE_OPERATORS_H
#define HLDO_SPHERE_OPERATORS_H

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
 * @brief Contains assembly classes for Hodge Laplacian operators, Dirac
 * Operators and the correspondig source problems
 */
namespace operators {}

}  // namespace projects::hldo_sphere

#endif  //  HLDO_SPHERE_OPERATORS_H
