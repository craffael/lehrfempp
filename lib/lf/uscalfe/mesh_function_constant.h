/**
 * @file
 * @brief Implements a \ref mesh_function "MeshFunction" which is constant on
 * the whole mesh.
 * @author Raffael Casagrande
 * @date   2018-12-16 07:29:52
 * @copyright MIT License
 */

#ifndef __d0f3b8f133da4af980ce21ffffdf719a
#define __d0f3b8f133da4af980ce21ffffdf719a
#include <vector>
#include "lf/mesh/mesh_interface.h"

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @ingroup mesh_function
 * @brief A \ref mesh_function "MeshFunction" which takes the same constant
 * value on the whole mesh.
 * @tparam R The type of the value.
 */
template <class R>
class MeshFunctionConstant {
 public:
  /**
   * @brief Create a new MeshFunctionConstant by passing the globally constant
   * value.
   * @param value The value that the MeshFunction should have globally.
   */
  explicit MeshFunctionConstant(R value) : value_(value) {}

  std::vector<R> operator()(const mesh::Entity& /*unused*/,
                            const Eigen::MatrixXd& local) const {
    return std::vector<R>(local.cols(), value_);
  }

 private:
  R value_;
};

}  // namespace lf::uscalfe

#endif  // __d0f3b8f133da4af980ce21ffffdf719a
