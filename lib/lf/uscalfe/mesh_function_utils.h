/**
 * @file
 * @brief Utility functions that act on mesh functions
 * @author Raffael Casagrande
 * @date   2019-01-13 09:04:29
 * @copyright MIT License
 */

#ifndef __faf48a519bd542b084b53cf30277a367
#define __faf48a519bd542b084b53cf30277a367

#include "mesh_function_traits.h"

namespace lf::uscalfe {

template <class MF, class QR_PROVIDER,
          class = std::enable_if_t<isMeshFunction<MF>>>
auto integrate(const lf::mesh::Mesh& mesh, const MF& mesh_function,
               QR_PROVIDER qr_provider, int codim = 0) {
  using return_t = MeshFunctionReturnType<MF>;
  return_t result;

  bool first = true;
  for (auto& e : mesh.Entities(codim)) {
    decltype(auto) qr = qr_provider(e);
    auto vals = mesh_function(e, qr.points());

    return_t temp;
    if constexpr (std::is_arithmetic_v<return_t>) {
      // Mesh function is scalar valued
      Eigen::Map<Eigen::Matrix<return_t, 1, Eigen::Dynamic>> vm(&vals[0], 1,
                                                                qr.NumPoints());
      temp = (vm * qr.Weights())(0, 0);
    } else if constexpr (std::is_convertible_v<return_t, Eigen::MatrixXd> ||
                         std::is_convertible_v<return_t, Eigen::MatrixXcd>) {
      // Mesh function is (Eigen) vector/matrix valued
      using Scalar = return_t::Scalar;

      if constexpr (return_t::RowsAtCompileTime == 1) {
        Eigen::Map<
            Eigen::Matrix<Scalar, return_t::SizeAtCompileTime, Eigen::Dynamic>>
            vm(&vals[0](0, 0), return_t::SizeAtCompileTime, vals.size());
        temp = vm * qr.Weights();
      } else if constexpr (return_t::ColsAtCompileTime == 1) {
        Eigen::Map<
            Eigen::Matrix<Scalar, return_t::SizeAtCompileTime, Eigen::Dynamic>>
            vm(&vals[0](0, 0), return_t::SizeAtCompileTime, vals.size());
        temp = (vm * qr.Weights()).transpose();
      } else {
        temp = qr.Weights()(0) * vals[0];
        for (int i = 1; i < qr.NumPoints(); ++i) {
          temp += qr.Weights()(i) * vals[i];
        }
      }
    } else {
      // Mesh function is otherwise valued -> don't use any optimizations.
      temp = qr.Weights()(0) * vals[0];
      for (int i = 1; i < qr.NumPoints(); ++i) {
        temp += qr.Weights()(i) * vals[i];
      }
    }

    if (first) {
      result = temp;
      first = false;
    } else {
      result += temp;
    }
  }

  return result;
}

}  // namespace lf::uscalfe

#endif  // __faf48a519bd542b084b53cf30277a367
