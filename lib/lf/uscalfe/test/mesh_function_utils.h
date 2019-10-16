/**
 * @file
 * @brief Some test utility functions related to mesh functions
 * @author Raffael Casagrande
 * @date   2019-01-13 08:08:40
 * @copyright MIT License
 */

#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe::test {

/** Checks whether two mesh functions are equal */
template <class A, class B>
void checkMeshFunctionEqual(const mesh::Mesh& m, A a, B b, int codim = 0) {
  using scalar_t = MeshFunctionReturnType<A>;
  static_assert(std::is_convertible_v<scalar_t, MeshFunctionReturnType<B>>);

  if constexpr (std::is_arithmetic_v<scalar_t>) {
    for (auto e : m.Entities(codim)) {
      auto ref_el = e->RefEl();
      auto qr = lf::quad::make_QuadRule(ref_el, 5);
      auto vals1 = a(*e, qr.Points());
      auto vals2 = b(*e, qr.Points());
      for (int i = 0; i < vals1.size(); ++i) {
        EXPECT_LT(vals1[i] - vals2[i], 1e-10) << "i=" << i;
      }
    }
  } else if constexpr (std::is_convertible_v<MeshFunctionReturnType<A>,
                                             Eigen::MatrixXd>) {
    for (auto e : m.Entities(codim)) {
      auto ref_el = e->RefEl();
      auto qr = lf::quad::make_QuadRule(ref_el, 5);
      auto vals1 = a(*e, qr.Points());
      auto vals2 = b(*e, qr.Points());
      for (int i = 0; i < vals1.size(); ++i) {
        EXPECT_LT((vals1[0] - vals2[0]).norm(), 1e-10);
      }
    }
  } else {
    for (auto e : m.Entities(codim)) {
      auto ref_el = e->RefEl();
      auto qr = lf::quad::make_QuadRule(ref_el, 5);
      auto vals1 = a(*e, qr.Points());
      auto vals2 = b(*e, qr.Points());
      for (int i = 0; i < vals1.size(); ++i) {
        EXPECT_EQ(vals1[i], vals2[i]) << "i=" << i;
      }
    }
  }
}

}  // namespace lf::uscalfe::test
