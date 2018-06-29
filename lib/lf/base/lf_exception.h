/**
 * @file
 * @brief Declares the LfException class.
 * @author Raffael Casagrande
 * @date   2018-06-29 05:58:49
 * @copyright MIT License
 */

#ifndef __d0efe6ebb86049268fc644633590ed83
#define __d0efe6ebb86049268fc644633590ed83
#include <exception>
#include <string>
#include <utility>

namespace lf::base {

/**
 * @brief A simple, general purpose exception class that is thrown by
 *        LehrFEM++ if something is wrong.
 */
class LfException : public std::exception {
 public:
  /**
   * @brief Create a new LfException with an error message.
   * @param what
   */
  LfException(std::string what) : what_(std::move(what)) {}

  char const* what() const noexcept override { return what_.c_str(); }

 private:
  std::string what_;
};

}  // namespace lf::base

#endif  // __d0efe6ebb86049268fc644633590ed83
