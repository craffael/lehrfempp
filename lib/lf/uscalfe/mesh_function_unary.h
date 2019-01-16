/**
 * @file
 * @brief overloads unary operators on mesh functions
 * @author Raffael Casagrande
 * @date   2019-01-13 07:32:39
 * @copyright MIT License
 */

#ifndef __b9b63bcccec548419a52fe0b06ffb3fc
#define __b9b63bcccec548419a52fe0b06ffb3fc

#include <lf/mesh/mesh.h>

namespace lf::uscalfe {

template <class OP, class MF>
class MeshFunctionUnary {
 public:
  MeshFunctionUnary(OP op, MF mf) : op_(std::move(op)), mf_(std::move(mf)) {}

  auto operator()(const mesh::Entity& e, const Eigen::MatrixXd& local) const {
    return op_(mf_(e, local), 0);
  }

 private:
  OP op_;
  MF mf_;
};

namespace internal {
struct UnaryOpMinus {
  // minus in front of a scalar type
  template <class U, class = std::enable_if_t<std::is_arithmetic_v<U>>>
  auto operator()(const std::vector<U>& u, int) const {
    Eigen::Map<const Eigen::Matrix<U, 1, Eigen::Dynamic>> um(&u[0], 1,
                                                             u.size());
    std::vector<U> result(u.size());
    Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>> rm(&result[0], 1, u.size());
    rm = -um;
    return result;
  }

  // minus in front of a fixed size matrix
  template <
      class S, int R, int C, int O, int MR, int MC,
      class = std::enable_if_t<R != Eigen::Dynamic && C != Eigen::Dynamic>>
  auto operator()(const std::vector<Eigen::Matrix<S, R, C, O, MR, MC>>& u,
                  int) const {
    if (R == 0 || C == 0) {
      // result vector is empty
      return u;
    }
    Eigen::Map<const Eigen::Matrix<S, 1, Eigen::Dynamic>> um(&u[0](0, 0), 1,
                                                             u.size() * R * C);
    std::vector<Eigen::Matrix<S, R, C, O, MR, MC>> result(u.size());
    Eigen::Map<Eigen::Matrix<S, 1, Eigen::Dynamic>> rm(&result[0](0, 0), 1,
                                                       u.size() * R * C);
    rm = -um;
    return result;
  }

  // minus in front of any object that supports the unary operator-
  template <class T>
  auto operator()(const std::vector<T>& u, long) const {
    std::vector<T> result(u.size());
    for (int i = 0; i < u.size(); ++i) {
      result[i] = -u[i];
    }
    return result;
  }
};

struct UnaryOpSquaredNorm {
  // squared norm of a scalar type
  template <class U, class = std::enable_if_t<std::is_arithmetic_v<U>>>
  auto operator()(const std::vector<U>& u, int) const {
    Eigen::Map<const Eigen::Matrix<U, 1, Eigen::Dynamic>> um(&u[0], 1,
                                                             u.size());
    std::vector<U> result(u.size());
    Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>> rm(&result[0], 1, u.size());
    rm = um.cwiseAbs2();
    return result;
  }

  // squared norm of a fixed size matrix
  template <
      class S, int R, int C, int O, int MR, int MC,
      class = std::enable_if_t<R != Eigen::Dynamic && C != Eigen::Dynamic>>
  auto operator()(const std::vector<Eigen::Matrix<S, R, C, O, MR, MC>>& u,
                  int) const {
    static_assert(R > 0 && C > 0,
                  "squaredNorm only supported for matrices with at least 1 row "
                  "and column");
    std::vector<double> result(u.size());
    if constexpr (C == 1) {
      Eigen::Map<const Eigen::Matrix<S, R, Eigen::Dynamic>> um(&u[0](0, 0), R,
                                                               u.size());
      Eigen::Map<Eigen::Matrix<S, 1, Eigen::Dynamic>> rm(&result[0], 1,
                                                         u.size());
      rm = um.cwiseAbs2().colwise().sum();
    } else if constexpr (R == 1) {
      Eigen::Map<const Eigen::Matrix<S, Eigen::Dynamic, C, Eigen::RowMajor>> um(
          &u[0](0, 0), u.size(), C);
      Eigen::Map<const Eigen::Matrix<S, Eigen::Dynamic, 1>> rm(&result[0],
                                                               u.size(), 1);
      rm = um.cwiseAbs2().rowwise().sum();
    } else {
      for (int i = 0; i < u.size(); ++i) {
        result[i] = u[i].squaredNorm();
      }
    }

    return result;
  }
};
}  // namespace internal

template <class A, class = std::enable_if_t<isMeshFunction<A>>>
auto operator-(const A& a) {
  return MeshFunctionUnary(internal::UnaryOpMinus{}, a);
}

template <class A, class = std::enable_if_t<isMeshFunction<A>>>
auto squaredNorm(const A& a) {
  return MeshFunctionUnary(internal::UnaryOpSquaredNorm{}, a);
}

}  // namespace lf::uscalfe

#endif  // __b9b63bcccec548419a52fe0b06ffb3fc
