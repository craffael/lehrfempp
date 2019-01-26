/**
 * @file
 * @brief Test the order of quadrature rules by testing it with monomials
 * @author Raffael Casagrande
 * @date   2018-08-19 06:54:02
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/quad/quad.h>
#include <boost/math/special_functions/factorials.hpp>

namespace lf::quad::test {

double integrate(QuadRule qr, std::vector<int> monomialCoefficients) {
  double result = 0;
  for (int i = 0; i < qr.Points().cols(); ++i) {
    double temp = 1;
    for (int j = 0; j < monomialCoefficients.size(); ++j) {
      temp *= std::pow(qr.Points()(j, i), monomialCoefficients[j]);
    }
    result += temp * qr.Weights()(i);
  }
  return result;
}

void checkQuadRule(QuadRule qr, double precision = 1e-12,
                   bool check_order_exact = true) {
  EXPECT_EQ(qr.Points().cols(), qr.Weights().size());
  EXPECT_EQ(qr.Points().rows(), qr.RefEl().Dimension());

  auto order = qr.Order();
  if (qr.RefEl() == base::RefEl::kSegment()) {
    for (int i = 0; i <= order; ++i) {
      // integrate x^i
      EXPECT_DOUBLE_EQ(integrate(qr, {i}), 1. / (1. + i))
          << "Failure for i = " << i;
    }
    // try integrate one order too high:
    EXPECT_GT(std::abs(integrate(qr, {static_cast<int>(order) + 1}) -
                       1. / (2. + order)),
              1e-10);
  } else if (qr.RefEl() == base::RefEl::kTria()) {
    // TRIA
    ///////////////////////////////////////////////////////////////////////////
    auto exact_value = [](int i, int j) {
      return boost::math::factorial<double>(i) *
             boost::math::factorial<double>(j + 1) /
             ((1 + j) * boost::math::factorial<double>(2 + i + j));
    };
    for (int i = 0; i <= order; ++i) {
      for (int j = 0; j <= order - i; ++j) {
        // integrate x^i y^j
        EXPECT_NEAR(integrate(qr, {i, j}) / exact_value(i, j), 1, precision);
      }
    }
    if (check_order_exact) {
      // Make sure that at least on of the order+1 polynomials is not integrated
      // correctly
      bool one_fails = false;
      for (int i = -1; i <= order; ++i) {
        if (std::abs(integrate(qr, {i + 1, static_cast<int>(order - i)}) -
                     exact_value(i + 1, order - i)) > 1e-12) {
          one_fails = true;
          break;
        }
      }
      EXPECT_TRUE(one_fails) << "order = " << (int)order;
    }

  } else if (qr.RefEl() == base::RefEl::kQuad()) {
    // QUAD
    ///////////////////////////////////////////////////////////////////////////
    for (int i = 0; i <= order; ++i) {
      for (int j = 0; j <= order; ++j) {
        // integrate x^i y^j
        EXPECT_DOUBLE_EQ(integrate(qr, {i, j}), 1. / ((1. + i) * (1. + j)));
      }
    }

    // make sure that not all of the higher polynomials integrate correctly:
    bool atLeastOneFails = false;
    for (int i = 0; i <= order + 1; ++i) {
      if (std::abs(integrate(qr, {static_cast<int>(order + 1), i}) -
                   1. / ((2. + order) * (1. + i))) > 1e-10) {
        atLeastOneFails = true;
        break;
      }
      if (std::abs(integrate(qr, {i, static_cast<int>(order + 1)}) -
                   1. / ((2. + order) * (1. + i))) > 1e-10) {
        atLeastOneFails = true;
        break;
      }
    }
    EXPECT_TRUE(atLeastOneFails);
  }
}

TEST(qr_IntegrationTest, Segment) {
  for (int i = 1; i < 10; ++i) {
    checkQuadRule(make_QuadRule(base::RefEl::kSegment(), i));
  }
}

TEST(qr_IntegrationTest, Quad) {
  checkQuadRule(make_QuadRule(base::RefEl::kQuad(), 1));
  checkQuadRule(make_QuadRule(base::RefEl::kQuad(), 2));
  checkQuadRule(make_QuadRule(base::RefEl::kQuad(), 3));
}

TEST(qr_IntegrationTest, Tria) {
  // make sure that also the tensor product versions are tested.
  for (int i = 1; i < 55; ++i) {
    checkQuadRule(make_QuadRule(base::RefEl::kTria(), i), 1e-12, i < 10);
  }
}

}  // namespace lf::quad::test
