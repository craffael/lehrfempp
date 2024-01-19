#ifndef HLDO_SPHERE_DEBUGGING_H
#define HLDO_SPHERE_DEBUGGING_H

namespace projects::hldo_sphere {
/**
 * @brief Contains debugging experiments mainly for the Dirac Opeartor
 *
 * This namespace was constructed to chase an error which turned out to be
 * a sign error in the file projects::hldo_sphere::operators::dirac_operator.cc
 *
 * Details about the experiments implemented in this section and the reasoning
 * behind, can be found in the
 * thesis `Hodge-Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * chapter 5.
 *
 */
namespace debugging {

/**
 * @brief Concatenate objects defining an operator<<(std::ostream&)
 * @param args A variadic pack of objects implementing
 * `operator<<(std::ostream&)`
 * @returns A string with the objects concatenated
 */
template <typename... Args>
static std::string concat(Args &&...args) {
  std::ostringstream ss;
  (ss << ... << args);
  return ss.str();
}

}  // namespace debugging

}  // namespace projects::hldo_sphere

#endif  //  HLDO_SPHERE_DEBUGGING_H
