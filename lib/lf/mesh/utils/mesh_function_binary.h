/**
 * @file
 * @brief Defines a couple of mesh functions for binary operations
 * @author Raffael Casagrande
 * @date   2019-01-13 04:33:05
 * @copyright MIT License
 */

#ifndef __9bad469d38e04e8ab67391ce50c2c480
#define __9bad469d38e04e8ab67391ce50c2c480
#include <Eigen/Eigen>
#include <type_traits>
#include <vector>
#include "mesh_function_traits.h"

namespace lf::mesh::utils {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @ingroup mesh_function
 * @brief A \ref mesh_function which combines two other \ref mesh_function
 * "mesh functions" using a binary operator (advanced use).
 * @tparam OP The type of operator that combines the mesh functions.
 * @tparam A The type of the lhs mesh function.
 * @tparam B The type of the rhs mesh function.
 *
 * ### Requirements for OP
 * The Operator `OP` must fulfill the following requirements:
 * - It must be moveable
 * - It should overload `operator()` as follows:
 * ```
 * template <class U, class V>
 * std::vector<Z> operator()(const std::vector<U>& u, const std::vector<V>& v,
 * int)
 * ```
 * where `U` is the MeshFunctionReturnType of the lhs MeshFunction, `V` is the
 * MeshFunctionReturnType of the rhs MeshFunction and `Z` is the type of the
 * mesh function `A OP B`.
 *
 * @note Usually there is no need to use MeshFunctionBinary directly. There are
 * a number of operator overloads which use MeshFunctionBinary internally.
 *
 * @note The last `int` argument of `operator()` can be used to prefer certain
 * overloads over others, e.g. `operator()(const std::vector<U>& u, const
 * std::vector<V>& int)` takes higher precedence than `operator()(const
 * std::vector<U>& u, const std::vector<V>& long)`
 *
 */
template <class OP, class A, class B>
class MeshFunctionBinary {
 public:
  /**
   * @brief Create a new MeshFunctionBinary
   * @param op The operator to apply
   * @param a The lhs mesh function
   * @param b The rhs mesh function.
   */
  MeshFunctionBinary(OP op, A a, B b)
      : op_(std::move(op)), a_(std::move(a)), b_(std::move(b)) {}

  /**
   * see \ref mesh_function for details.
   */
  auto operator()(const lf::mesh::Entity& e,
                  const Eigen::MatrixXd& local) const {
    return op_(a_(e, local), b_(e, local), 0);
  }

 private:
  OP op_;
  A a_;
  B b_;
};

/**
 * @brief Contains `OP` types (as used by MeshFunctionBinary) which are used by
 * the respective operator overloads (e.g. `operator+(...)`) to combine two mesh
 * functions.
 */
namespace internal {

/**
 * @brief Used together with MeshFunctionBinary (`OP` template type) to
 * represent the pointwise addition of two mesh functions.
 *
 * This struct contains multiple overloads of `operator()` which specialize for
 * certain cases, e.g. optimized implementations are provided for
 * `Eigen::Matrix` valued mesh functions. There is also a general overload of
 * `operator()` which works for any types that define `operator+`. But the more
 * specialized overloads are preferred over this one thanks to them having a 3rd
 * argument of type `int` whereas the general overload accepts a `long`.
 * (See also MeshFunctionBinary for an explanation)
 */
struct OperatorAddition {
  /**
   * @brief Addition of two scalar types (`std::is_arithmetic_v<...> == true`)
   *
   * @note The implementation of this method makes use of the fact that the
   * memory layout of the data of `std::vector<U>` and
   * `Eigen::Matrix<U,1,Eigen::Dynamic>` is the same. This allows us to
   * reinterpret the `std::vector<U>` as an `Eigen::Matrix` which in turn allows
   * us to use the Eigen optimized `operator+`
   */
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

  /**
   * @brief Addition of two `Eigen::Matrix` types.
   */
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
        LF_ASSERT_MSG(u[i].rows() == v[i].rows(),
                      "mismatch in #rows of matrices added to each other.");
        LF_ASSERT_MSG(u[i].cols() == v[i].cols(),
                      "mismatch in #cols of matrices added to each other.");
        result[i] = u[i] + v[i];
      }
      return result;
    } else {  // NOLINT
      // add two dynamic sized matrices:
      std::vector<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(u[i].rows() == v[i].rows(),
                      "mismatch in #rows of matrices added to each other.");
        LF_ASSERT_MSG(u[i].cols() == v[i].cols(),
                      "mismatch in #cols of matrices added to each other.");
        result.emplace_back(u[i] + v[i]);
      }
      return result;
    }
  }

  /**
   * @brief Addition of two `Eigen::Array` types.
   */
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Array<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Array<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) + S2(0));
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic &&
                  R2 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // add two static size eigen arrays to each other
      static_assert(R1 == R2, "cannot add arrays with different #rows.");
      static_assert(C1 == C2, "cannot add arrays with different #cols.");

      Eigen::Map<const Eigen::Array<S1, 1, Eigen::Dynamic>> um(
          &u[0](0, 0), 1, u.size() * R1 * C1);
      Eigen::Map<const Eigen::Array<S2, 1, Eigen::Dynamic>> vm(
          &v[0](0, 0), 1, v.size() * R1 * C1);

      std::vector<Eigen::Array<scalar_t, R1, C1>> result(u.size());
      Eigen::Map<Eigen::Array<scalar_t, 1, Eigen::Dynamic>> rm(
          &result[0](0, 0), 1, u.size() * R1 * C1);
      rm = um + vm;
      return result;
    }
    if constexpr ((R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) ||
                  (R2 != Eigen::Dynamic && C2 != Eigen::Dynamic)) {  // NOLINT
      // One of the arrays is fixed size:
      constexpr int fixedRows = std::max(R1, R2);
      constexpr int fixedCols = std::max(C1, C2);

      std::vector<Eigen::Array<scalar_t, fixedRows, fixedCols>> result(
          u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(u[i].rows() == v[i].rows(),
                      "mismatch in #rows of arrays added to each other.");
        LF_ASSERT_MSG(u[i].cols() == v[i].cols(),
                      "mismatch in #cols of arrays added to each other.");
        result[i] = u[i] + v[i];
      }
      return result;
    } else {  // NOLINT
      // add two dynamic sized arrays:
      std::vector<Eigen::Array<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(u[i].rows() == v[i].rows(),
                      "mismatch in #rows of arrays added to each other.");
        LF_ASSERT_MSG(u[i].cols() == v[i].cols(),
                      "mismatch in #cols of arrays added to each other.");
        result.emplace_back(u[i] + v[i]);
      }
      return result;
    }
  }

  /**
   * @brief addition of arbitrary types supporting  `operator+`:
   */
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

/**
 * @brief Used together with MeshFunctionBinary (`OP` template type) to
 * represent the pointwise subtraction of two mesh functions.
 *
 * This struct contains multiple overloads of `operator()` which specialize for
 * certain cases, e.g. optimized implementations are provided for
 * `Eigen::Matrix` valued mesh functions. There is also a general overload of
 * `operator()` which works for any types that define `operator+`. But the more
 * specialized overloads are preferred over this one thanks to them having a 3rd
 * argument of type `int` whereas the general overload accepts a `long`.
 * (See also MeshFunctionBinary for an explanation)
 */
struct OperatorSubtraction {
  /**
   * @brief Subtraction of two scalar types (`std::is_arithmetic_v<...> ==
   * true`)
   *
   * @note The implementation of this method makes use of the fact that the
   * memory layout of the data of `std::vector<U>` and
   * `Eigen::Matrix<U,1,Eigen::Dynamic>` is the same. This allows us to
   * reinterpret the `std::vector<U>` as an `Eigen::Matrix` which in turn allows
   * us to use the Eigen optimized `operator-`
   */
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

  /**
   * @brief Subtraction of two `Eigen::Matrix` types.
   */
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Matrix<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Matrix<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) + S2(0));
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic &&
                  R2 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // subtract two static size eigen matrices from each other
      static_assert(R1 == R2, "cannot subtract matrices with different #rows.");
      static_assert(C1 == C2, "cannot subtract matrices with different #cols.");

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
      for (std::size_t i = 0; i < u.size(); ++i) {
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

  /**
   * @brief Subtraction of two `Eigen::Array` types.
   */
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Array<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Array<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) - S2(0));
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic &&
                  R2 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // subtract two static size eigen arrays to each other
      static_assert(R1 == R2, "cannot subtract arrays with different #rows.");
      static_assert(C1 == C2, "cannot subtract arrays with different #cols.");

      Eigen::Map<const Eigen::Array<S1, 1, Eigen::Dynamic>> um(
          &u[0](0, 0), 1, u.size() * R1 * C1);
      Eigen::Map<const Eigen::Array<S2, 1, Eigen::Dynamic>> vm(
          &v[0](0, 0), 1, v.size() * R1 * C1);

      std::vector<Eigen::Array<scalar_t, R1, C1>> result(u.size());
      Eigen::Map<Eigen::Array<scalar_t, 1, Eigen::Dynamic>> rm(
          &result[0](0, 0), 1, u.size() * R1 * C1);
      rm = um - vm;
      return result;
    }
    if constexpr ((R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) ||
                  (R2 != Eigen::Dynamic && C2 != Eigen::Dynamic)) {  // NOLINT
      // One of the arrays is fixed size:
      constexpr int fixedRows = std::max(R1, R2);
      constexpr int fixedCols = std::max(C1, C2);

      std::vector<Eigen::Array<scalar_t, fixedRows, fixedCols>> result(
          u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of arrays subtracted from each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of arrays subtracted from each other.");
        result[i] = u[i] - v[i];
      }
      return result;
    } else {  // NOLINT
      // subtract two dynamic sized arrays:
      std::vector<Eigen::Array<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of arrays subtracted from each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of arrays subtracted from each other.");
        result.emplace_back(u[i] - v[i]);
      }
      return result;
    }
  }

  /**
   *@brief subtraction of arbitrary types supporting - operator:
   */
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

/**
 * @brief Used together with MeshFunctionBinary (`OP` template type) to
 * represent the pointwise product of two mesh functions.
 *
 * This struct contains multiple overloads of `operator()` which specialize for
 * certain cases, e.g. optimized implementations are provided for
 * `Eigen::Matrix` valued mesh functions. There is also a general overload of
 * `operator()` which works for any types that define `operator+`. But the more
 * specialized overloads are preferred over this one thanks to them having a 3rd
 * argument of type `int` whereas the general overload accepts a `long`.
 * (See also MeshFunctionBinary for an explanation)
 */
struct OperatorMultiplication {
  /**
   * @brief Multiplication of two scalar types (`std::is_arithmetic_v<...> ==
   * true`)
   *
   * @note The implementation of this method makes use of the fact that the
   * memory layout of the data of `std::vector<U>` and
   * `Eigen::Array<U,1,Eigen::Dynamic>` is the same. This allows us to
   * reinterpret the `std::vector<U>` as an `Eigen::Array` which in turn allows
   * us to use the Eigen optimized `operator*`
   */
  template <class U, class V,
            class = typename std::enable_if<std::is_arithmetic_v<U> &&
                                            std::is_arithmetic_v<V>>::type>
  auto operator()(const std::vector<U>& u, const std::vector<V>& v, int
                  /*unused*/) const {
    Eigen::Map<const Eigen::Array<U, 1, Eigen::Dynamic>> um(&u[0], 1, u.size());
    Eigen::Map<const Eigen::Array<U, 1, Eigen::Dynamic>> vm(&v[0], 1, v.size());
    std::vector<decltype(U(0) * V(0))> result(u.size());
    Eigen::Map<Eigen::Array<decltype(U(0) * V(0)), 1, Eigen::Dynamic>> rm(
        &result[0], 1, u.size());
    rm = um * vm;
    return result;
  }

  // multiplication of fixed size eigen matrices
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Matrix<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Matrix<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) * S2(0));
    if constexpr (R1 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // The result is fixed size
      static_assert(C1 == Eigen::Dynamic || R2 == Eigen::Dynamic || C1 == R2,
                    "#cols of lhs doesn't match #rows of rhs");

      std::vector<Eigen::Matrix<scalar_t, R1, C2>> result(u.size());
      for (std::size_t i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(u[i].cols() == v[i].rows(),
                      "#cols of lhs doesn't match #rows of rhs");
        result[i] = u[i] * v[i];
      }
      return result;
    } else {  // NOLINT
      // multiply dynamic sized matrices
      std::vector<Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(u[i].cols() == v[i].rows(),
                      "#cols of lhs doesn't match #rows of rhs");
        result.emplace_back(u[i] * v[i]);
      }
      return result;
    }
  }

  /**
   * @brief Multiplication of two `Eigen::Array` types.
   */
  template <class S1, int R1, int C1, int O1, int MR1, int MC1, class S2,
            int R2, int C2, int O2, int MR2, int MC2>
  auto operator()(const std::vector<Eigen::Array<S1, R1, C1, O1, MR1, MC1>>& u,
                  const std::vector<Eigen::Array<S2, R2, C2, O2, MR2, MC2>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(S1(0) * S2(0));
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic &&
                  R2 != Eigen::Dynamic && C2 != Eigen::Dynamic) {  // NOLINT
      // multiply two static size eigen arrays to each other
      static_assert(R1 == R2, "cannot multiply arrays with different #rows.");
      static_assert(C1 == C2, "cannot multiply arrays with different #cols.");

      Eigen::Map<const Eigen::Array<S1, 1, Eigen::Dynamic>> um(
          &u[0](0, 0), 1, u.size() * R1 * C1);
      Eigen::Map<const Eigen::Array<S2, 1, Eigen::Dynamic>> vm(
          &v[0](0, 0), 1, v.size() * R1 * C1);

      std::vector<Eigen::Array<scalar_t, R1, C1>> result(u.size());
      Eigen::Map<Eigen::Array<scalar_t, 1, Eigen::Dynamic>> rm(
          &result[0](0, 0), 1, u.size() * R1 * C1);
      rm = um * vm;
      return result;
    }
    if constexpr ((R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) ||
                  (R2 != Eigen::Dynamic && C2 != Eigen::Dynamic)) {  // NOLINT
      // One of the arrays is fixed size:
      constexpr int fixedRows = std::max(R1, R2);
      constexpr int fixedCols = std::max(C1, C2);

      std::vector<Eigen::Array<scalar_t, fixedRows, fixedCols>> result(
          u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of arrays multiplied with each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of arrays multiplied with each other.");
        result[i] = u[i] * v[i];
      }
      return result;
    } else {  // NOLINT
      // multiply two dynamic sized arrays:
      std::vector<Eigen::Array<scalar_t, Eigen::Dynamic, Eigen::Dynamic>>
          result;
      result.reserve(u.size());
      for (int i = 0; i < u.size(); ++i) {
        LF_ASSERT_MSG(
            u[i].rows() == v[i].rows(),
            "mismatch in #rows of arrays multiplied with each other.");
        LF_ASSERT_MSG(
            u[i].cols() == v[i].cols(),
            "mismatch in #cols of arrays multiplied with each other.");
        result.emplace_back(u[i] * v[i]);
      }
      return result;
    }
  }

  // multiplication of a scalar with matrix
  template <class U, class S1, int R1, int C1, int O1, int MR1, int MC1,
            class = std::enable_if_t<std::is_arithmetic_v<U>>>
  auto operator()(const std::vector<U>& u,
                  const std::vector<Eigen::Matrix<S1, R1, C1, O1, MR1, MC1>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(u[0] * v[0](0, 0));
    std::vector<Eigen::Matrix<scalar_t, R1, C1>> result(u.size());
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) {
      // result is a static sized matrix:
      Eigen::Map<const Eigen::Array<S1, R1 * C1, Eigen::Dynamic>> vm(
          v[0].data(), R1 * C1, v.size());
      Eigen::Map<Eigen::Array<scalar_t, R1 * C1, Eigen::Dynamic>> rm(
          result[0].data(), R1 * C1, v.size());
      Eigen::Map<const Eigen::Array<U, 1, Eigen::Dynamic>> um(u.data(), 1,
                                                              u.size());

      rm = vm * um.template replicate<R1 * C1, 1>();
    } else {
      // result is not static sized:
      for (std::size_t i = 0; i < u.size(); ++i) {
        result[i] = u[i] * v[i];
      }
    }
    return result;
  }

  // multiplication of matrix with scalar (other way around)
  template <class U, class S1, int R1, int C1, int O1, int MR1, int MC1,
            class = std::enable_if_t<std::is_arithmetic_v<U>>>
  auto operator()(const std::vector<Eigen::Matrix<S1, R1, C1, O1, MR1, MC1>>& v,
                  const std::vector<U>& u, int /*unused*/) const {
    return operator()(u, v, 0);
  }

  // multiplication of a scalar with array
  template <class U, class S1, int R1, int C1, int O1, int MR1, int MC1,
            class = std::enable_if_t<std::is_arithmetic_v<U>>>
  auto operator()(const std::vector<U>& u,
                  const std::vector<Eigen::Array<S1, R1, C1, O1, MR1, MC1>>& v,
                  int /*unused*/) const {
    using scalar_t = decltype(u[0] * v[0](0, 0));
    std::vector<Eigen::Array<scalar_t, R1, C1>> result(u.size());
    if constexpr (R1 != Eigen::Dynamic && C1 != Eigen::Dynamic) {
      // result is a static sized array:
      Eigen::Map<const Eigen::Array<S1, R1 * C1, Eigen::Dynamic>> vm(
          v[0].data(), R1 * C1, v.size());
      Eigen::Map<Eigen::Array<scalar_t, R1 * C1, Eigen::Dynamic>> rm(
          result[0].data(), R1 * C1, v.size());
      Eigen::Map<const Eigen::Array<U, 1, Eigen::Dynamic>> um(u.data(), 1,
                                                              u.size());

      rm = vm * um.template replicate<R1 * C1, 1>();
    } else {
      // result is not static sized:
      for (std::size_t i = 0; i < u.size(); ++i) {
        result[i] = u[i] * v[i];
      }
    }
    return result;
  }

  // multiplication of array with scalar (other way around)
  template <class U, class S1, int R1, int C1, int O1, int MR1, int MC1,
            class = std::enable_if_t<std::is_arithmetic_v<U>>>
  auto operator()(const std::vector<Eigen::Array<S1, R1, C1, O1, MR1, MC1>>& v,
                  const std::vector<U>& u, int /*unused*/) const {
    return operator()(u, v, 0);
  }

  // multiplication of arbitrary types supporting operator*:
  template <class U, class V>
  auto operator()(const std::vector<U>& u, const std::vector<V>& v,
                  long /*unused*/) const {
    std::vector<decltype(u[0] * v[0])> result;
    result.reserve(u.size());
    for (int i = 0; i < result.size(); ++i) {
      result.emplace_back(u[i] * v[i]);
    }
    return result;
  }
};

}  // namespace internal

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Add's two \ref mesh_function "mesh functions"
 * @relatesalso lf::uscalfe::MeshFunctionBinary
 * @tparam A Type of the lhs \ref mesh_function
 * @tparam B Type of the rhs \ref mesh_function
 * @param a the lhs \ref mesh_function
 * @param b the rhs \ref mesh_function
 * @return `a + b`, i.e. a \ref mesh_function "new mesh function" which
 * represents the pointwise addition of `a` and `b`
 *
 * @note the two \ref mesh_function "mesh functions" `a` and `b` should produce
 * the same type of values, e.g. both should be scalar valued or both
 * matrix/vector/array valued.
 *
 * @note If you want to apply this operator overload to \ref mesh_function "mesh
 * functions" which do not reside in the `lf::uscalfe` namespace, it will not be
 * found by Argument Dependent Lookup. You can get around this by explictly
 * importing the operator overload: `using lf::uscalfe::operator+;`
 *
 * ### Example
 * @snippet mesh_function_binary.cc one_trig
 */
template <class A, class B,
          class = std::enable_if_t<isMeshFunction<A> && isMeshFunction<B>>>
auto operator+(const A& a, const B& b) {
  return MeshFunctionBinary(internal::OperatorAddition{}, a, b);
}

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Subtracts two \ref mesh_function "mesh functions"
 * @relatesalso lf::uscalfe::MeshFunctionBinary
 * @tparam A Type of the lhs \ref mesh_function
 * @tparam B Type of the rhs \ref mesh_function
 * @param a the lhs \ref mesh_function
 * @param b the rhs \ref mesh_function
 * @return `a - b`, i.e. a \ref mesh_function "new mesh function" which
 * represents the pointwise difference of `a` minus `b`
 *
 * @note the two \ref mesh_function "mesh functions" `a` and `b` should produce
 * the same type of values, e.g. both should be scalar valued or both
 * matrix/vector valued.
 *
 * @note If you want to apply this operator overload to \ref mesh_function "mesh
 * functions" which do not reside in the `lf::uscalfe` namespace, it will not be
 * found by Argument Dependent Lookup (ADL). You can get around this by
 * explicitly importing the operator overload: `using lf::uscalfe::operator-;`
 *
 * ### Example
 * @snippet mesh_function_binary.cc subtract
 */
template <class A, class B,
          class = std::enable_if_t<isMeshFunction<A> && isMeshFunction<B>>>
auto operator-(const A& a, const B& b) {
  return MeshFunctionBinary(internal::OperatorSubtraction{}, a, b);
}

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @relatesalso lf::uscalfe::MeshFunctionBinary
 * @brief Multiply two \ref mesh_function "mesh functions" with each other.
 * @tparam A The type of the lhs \ref mesh_function
 * @tparam B The type of the rhs \ref mesh_function
 * @param a the lhs \ref mesh_function
 * @param b the rhs \ref mesh_function
 * @return Return `a*b`, i.e. a \ref mesh_function "new mesh function" which
 * represents the pointwise product of `a` and `b`.
 *
 * @note if `A` and `B` are `Eigen::Matrix` valued, the resulting \ref
 * mesh_function will also be Matrix/Vector valued and will represent the Matrix
 * product of `a` and `b`.
 *
 * @note If you want to apply this operator overload to \ref mesh_function "mesh
 * functions" which do not reside in the `lf::uscalfe` namespace, it will not be
 * found by Argument Dependent Lookup (ADL). You can get around this by
 * explicitly importing the operator overload: `using lf::uscalfe::operator*;`
 *
 * ### Example
 * @snippet mesh_function_binary.cc product
 */
template <class A, class B,
          class = std::enable_if_t<isMeshFunction<A> && isMeshFunction<B>>>
auto operator*(const A& a, const B& b) {
  return MeshFunctionBinary(internal::OperatorMultiplication{}, a, b);
}

}  // namespace lf::mesh::utils

#endif  // __9bad469d38e04e8ab67391ce50c2c480
