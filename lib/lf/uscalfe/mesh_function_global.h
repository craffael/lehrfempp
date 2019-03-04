/**
 * @file
 * @brief Transforms an function expressed in global coordinates into a
 * MeshFunction.
 * @author Raffael Casagrande
 * @date   2018-12-15 09:29:05
 * @copyright MIT License
 */

#ifndef __1323b89d2f584608940fcb37657142cc
#define __1323b89d2f584608940fcb37657142cc
#include <vector>

#include <lf/mesh/mesh.h>

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @ingroup mesh_function
 * @brief MeshFunction wrapper for a simple function of physical coordinates
 *
 * @tparam F a functor type, offering an evaluation operator that acccepts
 *           a single coordinate vector.
 *
 * MeshFunctionGlobal essentially wraps a [function
 * object](https://en.wikipedia.org/wiki/Function_object) which represents a
 * function defined on the mesh that depends only on the global coordinates.
 * An example is e.g. the function \f$ \vec{x} \mapsto \norm{x}^2 \f$.
 *
 * ### Requirements for `F`
 * `F` is a function object that should overload the call operator as follows:
 * ```
 * auto operator()(const Eigen::Vector2d& x) const
 * ```
 * which should return the value of the mesh function at the global coordinates
 * `x`. The return type of the call operator can in principle be anything, but
 * usually it is one of:
 * - `double` for a scalar valued mesh function
 * - `std::complex<double>` for complex valued mesh function
 * - `Eigen::Vector2d` for a tensor valued mesh function.
 *
 * @note For `DimGlobal==3`, the call operator should of course accept a
 * `Eigen::Vector3d`.
 */
template <class F>
class MeshFunctionGlobal {
  using F_return_type =
      decltype(std::declval<F>()(std::declval<Eigen::Vector2d>()));

 public:
  MeshFunctionGlobal(const MeshFunctionGlobal&) = default;
  MeshFunctionGlobal(MeshFunctionGlobal&&) noexcept = default;
  MeshFunctionGlobal& operator=(const MeshFunctionGlobal&) = delete;
  MeshFunctionGlobal& operator=(MeshFunctionGlobal&&) = delete;

  explicit MeshFunctionGlobal(F f) : f_(std::move(f)) {}

  /**
   * @brief MeshFunction compliant evaluation operator
   */
  std::vector<F_return_type> operator()(const mesh::Entity& e,
                                        const Eigen::MatrixXd& local) const {
    LF_ASSERT_MSG(e.RefEl().Dimension() == local.rows(),
                  "mismatch between entity dimension and local.rows()");
    std::vector<F_return_type> result;
    result.reserve(local.cols());
    auto global_points = e.Geometry()->Global(local);
    for (base::dim_t i = 0; i < local.cols(); ++i) {
      result.push_back(f_(global_points.col(i)));
    }
    return result;
  }

  virtual ~MeshFunctionGlobal() = default;

 private:
  F f_;
};

}  // namespace lf::uscalfe

#endif  // __1323b89d2f584608940fcb37657142cc
