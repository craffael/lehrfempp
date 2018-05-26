#ifndef __b58a30690f684a2b9422606fd526d0f4
#define __b58a30690f684a2b9422606fd526d0f4

#include <exception>
#include <utility>
#include <string>

namespace lf::base {

/**
 * @brief Thrown to signal that an argument passed to a function had the
 * wrong (polymorphic) type.
 */
class InvalidTypeException : public std::exception {
 private:
  std::string what_;

 public:
  InvalidTypeException() = default;

  explicit InvalidTypeException(std::string message)
      : what_(std::move(message)) {}

  InvalidTypeException(const InvalidTypeException& other) = default;
};

}  // namespace lf::base

#endif  // __b58a30690f684a2b9422606fd526d0f4
