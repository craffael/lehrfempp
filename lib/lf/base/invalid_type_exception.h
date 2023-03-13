#ifndef INCGb58a30690f684a2b9422606fd526d0f4
#define INCGb58a30690f684a2b9422606fd526d0f4

#include <exception>
#include <string>
#include <utility>

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
  InvalidTypeException(const InvalidTypeException&) = default;
  InvalidTypeException(InvalidTypeException&&) = default;

  explicit InvalidTypeException(std::string message)
      : what_(std::move(message)) {}

  InvalidTypeException& operator=(const InvalidTypeException&) = delete;
  InvalidTypeException& operator=(InvalidTypeException&&) = delete;

  ~InvalidTypeException() override = default;
};

}  // namespace lf::base

#endif  // INCGb58a30690f684a2b9422606fd526d0f4
