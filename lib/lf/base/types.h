/**
 * @file
 * @brief Defines common LehrFEM++ types
 * @author Raffael Casagrande
 * @date   2024-01-05 16:34
 * @copyright MIT License
 */

#ifndef INCG_e1ab6004cad5460f9baa3153c65f157a
#define INCG_e1ab6004cad5460f9baa3153c65f157a

namespace lf::base {
/** @defgroup lftypes Common Typedefs
 * @brief various integral types meant to enhance readability of the code
 */
/** @{ */
/**
 * @brief general type for variables related to size of arrays
 */
using size_type = unsigned int;
/**
 * @brief type for global index of mesh entities (nodes, edges, cells)
 */
using glb_idx_t = unsigned int;
/**
 * @brief type for local indices of sub-entities
 */
using sub_idx_t = unsigned int;
/**
 * @brief type for dimensions and co-dimensions and numbers derived from them
 */
using dim_t = unsigned int;
/**
 * @brief Index flagged as invalid
 */
const unsigned int kIdxNil = static_cast<unsigned int>(-1);
/** @} */

constexpr double kPi = 3.14159265358979323846;
}  // namespace lf::base

#endif  // INCG_e1ab6004cad5460f9baa3153c65f157a