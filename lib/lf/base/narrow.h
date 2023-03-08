/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Helper function to make casts explicit
 * @author Raffael Casagrande
 * @date March 2023
 * @copyright MIT License
 */

#ifndef INCG_1944399668484015bcf9f347522a68c3
#define INCG_1944399668484015bcf9f347522a68c3
#include <type_traits>

namespace lf::base {

/**
 * @brief cast `from` to `TO` using `static_cast`. An assert will fail if this
 * would change the value (e.g. because of overflow when an `std::int64_t` is
 * converted to a `std::int32_t`)
 * @tparam TO Type to which we want to convert
 * @tparam FROM Type from which we convert
 * @param from original value
 * @return `static_cast<TO>(std::forward<FROM>(from))`
 */
template <class TO, class FROM>
auto narrow(FROM from) noexcept -> TO {
  static_assert(std::is_arithmetic_v<TO>, "TO must be an arithmetic type");
  const TO to = static_cast<TO>(std::forward<FROM>(from));
  if constexpr (std::is_arithmetic_v<TO>) {
    LF_ASSERT_MSG(static_cast<FROM>(to) == from &&
                      (std::is_signed_v<FROM> == std::is_signed_v<TO> ||
                       ((to < TO{}) == (from < FROM{}))),
                  "Narrowing error.");
  } else {
    LF_ASSERT_MSG(static_cast<FROM>(to) == from, "Narrowing error.");
  }
  return to;
}

}  // namespace lf::base

#endif  // INCG_1944399668484015bcf9f347522a68c3