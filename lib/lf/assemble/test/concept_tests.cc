/**
 * @file
 * @brief Test concepts defined in assemble/concepts.h
 * @author Raffael Casagrande
 * @date   12.1.2024
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>

namespace lf::assemble::tests {

// Make sure the Archetypes really fulfills the concept.
static_assert(EntityMatrixProvider<EntityMatrixProviderAT<double>>);
static_assert(
    EntityMatrixProvider<EntityMatrixProviderAT<std::complex<double>>>);

static_assert(EntityVectorProvider<EntityVectorProviderAT<double>>);
static_assert(
    EntityVectorProvider<EntityVectorProviderAT<std::complex<double>>>);

}  // namespace lf::assemble::tests
