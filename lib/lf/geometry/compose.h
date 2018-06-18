/**
 * @file
 * @brief Contains the lf::geometry::Compose() global function.
 * @author Raffael Casagrande
 * @date   2018-06-17 06:24:51
 * @copyright MIT License
 */

#ifndef __ae8e3269a3cd49398e40713cd251f128
#define __ae8e3269a3cd49398e40713cd251f128
#include "geometry_interface.h"

namespace lf::geometry {

/**
 * @brief Create a new geometry object that is the composition of
 *        \f$ \Phi_a \circ \Phi_b \f$
 * @param a The mapping \f$ \Phi_a \f$
 * @param b The mapping \f$ \Phi_b \f$
 * @return The composition of \f$\Phi_a \circ \Phi_b\f$
 */
std::unique_ptr<Geometry> Compose(std::unique_ptr<Geometry>&& a,
                                  std::unique_ptr<Geometry>&& b);

}  // namespace lf::geometry

#endif  // __ae8e3269a3cd49398e40713cd251f128
