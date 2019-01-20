/**
 * @file
 * @brief Defines a couple of mesh functions for binary operations
 * @author Raffael Casagrande
 * @date   2019-01-13 04:33:05
 * @copyright MIT License
 */

#ifndef __9bad469d38e04e8ab67391ce50c2c480
#define __9bad469d38e04e8ab67391ce50c2c480

#include "lf/mesh/mesh_interface.h"

namespace lf::uscalfe {

/**
 * @brief A mesh function which combines two other mesh functions using a binary
 * operator.
 * @tparam OP The type of operator that combines the mesh functions.
 * @tparam A The type of the lhs mesh function.
 * @tparam B The type of the rhs mesh function.
 *
 * # Requirements for OP
 *
 */
template <class OP, class A, class B>
class MeshFunctionBinary {
 public:
  MeshFunctionBinary(OP op, A a, B b)
      : op_(std::move(op)), a_(std::move(a)), b_(std::move(b)) {}

  auto operator()(const lf::mesh::Entity& e,
                  const Eigen::MatrixXd& local) const {
    return op_(a_(e, local), b_(e, local), 0);
  }

 private:
  OP op_;
  A a_;
  B b_;
};

namespace internal {
struct OperatorAddition {
  // addition of scalar types:
  template <class U, class V,
            class = typename std::enable_if<std::is_arithmetic_v<U> &&
                                            std::is_arithmetic_v<V>>::type>
  auto operator()(const std::vector<U>& u, const std::vector<V>& v, int
                  /*unused*/) const {
    Eigen::Map<const Eigen::Matrix<U, 1, Eigen::Dynamic>> um(&u[0], 1,
                                                             u.size());
    Eigen::Map<const Eigen::Matrix<U, 1, Eigen::Dynamic>> vm(&v[0], 1,
                                                             v.size());
    std::vector<decltype(U(0) + V(0))> result(u.size());
    Eigen::Map<Eigen::Matrix<decltype(U(0) + V(0)), 1, Eigen::Dynamic>> rm(
        &result[0], 1, u.size());
    rm = um + vm;
    return result;
  }

  // Addition of static sized eigen matrices
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Matrix<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Matrix<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) + S2(0));
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic &&
                  R2 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // add two static size eigen matrices to each other
      static_assert(R1 == R2, "cannot add matrices with different #rows.");
      static_assert(C1 == C2, "cannot add matrices with different #cols.");

      Eigen::Map<const Eigen::Matrix<S1, 1, Eigen::Dynamic>> um(
          &u[0](0, 0), 1, u.size() * R1 * C1);
      Eigen::Map<const Eigen::Matrix<S2, 1, Eigen::Dynamic>> vm(
          &v[0](0, 0), 1, v.size() * R1 * C1);

      std::vector<Eigen::Matrix<scalar_t, R1, C1>> result(u.size());
      Eigen::Map<Eigen::Matrix<scalar_t, 1, Eigen::Dynamic>> rm(
          &result[0](0, 0), 1, u.size() * R1 * C1);
      rm = um + vm;
      return result;
    }
    if constexpr ((R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) ||
                  (R2 != Eigen::Dynamic && C2 != Eigen::Dynamic)) {  // NOLINT
      // One of the matrices is fixed size:
      constexpr int fixedRows = std::max(R1, R2);
      constexpr int fixedCols = std::max(C1, C2);

      std::vector<Eigen::Matrix<scalar_t, fixedRows, fixedCols>> result(
          u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of matrices subtracted from each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of matrices subtracted from each other.");
        result[i] = u[i] + v[i];
      }
      return result;
    } else {  // NOLINT
      // subtract two dynamic sized matrices from each other:
      std::vector<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of matrices subtracted from each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of matrices subtracted from each other.");
        result.emplace_back(u[i] + v[i]);
      }
      return result;
    }
  }

  // addition of arbitrary types supporting + operator:
  template <class U, class V>
  auto operator()(const std::vector<U>& u, const std::vector<V>& v,
                  long /*unused*/) const {
    std::vector<decltype(u[0] + v[0])> result(u.size());
    for (int i = 0; i < result.size(); ++i) {
      result[i] = u[i] + v[i];
    }
    return result;
  }
};

struct OperatorSubtraction {
  // subtraction of scalar types:
  template <class U, class V,
            class = typename std::enable_if<std::is_arithmetic_v<U> &&
                                            std::is_arithmetic_v<V>>::type>
  auto operator()(const std::vector<U>& u, const std::vector<V>& v, int
                  /*unused*/) const {
    Eigen::Map<const Eigen::Matrix<U, 1, Eigen::Dynamic>> um(&u[0], 1,
                                                             u.size());
    Eigen::Map<const Eigen::Matrix<U, 1, Eigen::Dynamic>> vm(&v[0], 1,
                                                             v.size());
    std::vector<decltype(U(0) + V(0))> result(u.size());
    Eigen::Map<Eigen::Matrix<decltype(U(0) + V(0)), 1, Eigen::Dynamic>> rm(
        &result[0], 1, u.size());
    rm = um - vm;
    return result;
  }

  // subtraction of fixed size eigen matrices
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Matrix<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Matrix<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) + S2(0));
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic &&
                  R2 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // subtract two static size eigen matrices from each other
      static_assert(R1 == R2, "cannot add matrices with different #rows.");
      static_assert(C1 == C2, "cannot add matrices with different #cols.");

      Eigen::Map<const Eigen::Matrix<S1, 1, Eigen::Dynamic>> um(
          &u[0](0, 0), 1, u.size() * R1 * C1);
      Eigen::Map<const Eigen::Matrix<S2, 1, Eigen::Dynamic>> vm(
          &v[0](0, 0), 1, v.size() * R1 * C1);

      std::vector<Eigen::Matrix<scalar_t, R1, C1>> result(u.size());
      Eigen::Map<Eigen::Matrix<scalar_t, 1, Eigen::Dynamic>> rm(
          &result[0](0, 0), 1, u.size() * R1 * C1);
      rm = um - vm;
      return result;
    }
    if constexpr ((R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) ||
                  (R2 != Eigen::Dynamic && C2 != Eigen::Dynamic)) {  // NOLINT
      // One of the matrices is fixed size:
      constexpr int fixedRows = std::max(R1, R2);
      constexpr int fixedCols = std::max(C1, C2);

      std::vector<Eigen::Matrix<scalar_t, fixedRows, fixedCols>> result(
          u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of matrices subtracted from each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of matrices subtracted from each other.");
        result[i] = u[i] - v[i];
      }
      return result;
    } else {  // NOLINT
      // subtract two dynamic sized matrices from each other:
      std::vector<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of matrices subtracted from each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of matrices subtracted from each other.");
        result.emplace_back(u[i] - v[i]);
      }
      return result;
    }
  }

  // subtraction of arbitrary types supporting - operator:
  template <class U, class V>
  auto operator()(const std::vector<U>& u, const std::vector<V>& v,
                  long /*unused*/) const {
    std::vector<decltype(u[0] + v[0])> result(u.size(),
                                              decltype(u[0] + v[0])());
    for (int i = 0; i < result.size(); ++i) {
      result[i] = u[i] - v[i];
    }
    return result;
  }
};

}  // namespace internal

template <class A, class B,
          class = std::enable_if_t<isMeshFunction<A> && isMeshFunction<B>>>
auto operator+(const A& a, const B& b) {
  return MeshFunctionBinary(internal::OperatorAddition{}, a, b);
}

template <class A, class B,
          class = std::enable_if_t<isMeshFunction<A> && isMeshFunction<B>>>
auto operator-(const A& a, const B& b) {
  return MeshFunctionBinary(internal::OperatorSubtraction{}, a, b);
}

}  // namespace lf::uscalfe

#endif  // __9bad469d38e04e8ab67391ce50c2c480
