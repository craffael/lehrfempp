/** @file runtime_test.cc
 *  @brief Some tests for runtime behavior of certain C++ constructs
 */

#include <iostream>

#include "lf/base/base.h"

static const int N = 10;

static double stat_tmp[N];  // NOLINT

std::vector<double> getData_VEC(double offset) {
  std::vector<double> tmp(10);
  for (int j = 0; j < N; j++) {
    tmp[j] = 1.0 / (j + offset);
  }
  return tmp;
}

// NOLINTNEXTLINE
void getData_REF(double offset, std::vector<double> &res) {
  for (int j = 0; j < N; j++) {
    res[j] = 1.0 / (j + offset);
  }
}

nonstd::span<const double> getData_SPAN(double offset) {
  for (int j = 0; j < N; j++) {
    stat_tmp[j] = 1.0 / (j + offset);
  }
  return {static_cast<double *>(stat_tmp), (stat_tmp + N)};
}

int main(int /*argc*/, const char * /*unused*/[]) {
  std::cout << "Runtime test for range access" << '\n';

  const long int reps = 100000000L;

  std::cout << "I. Returning std::vector's" << '\n';
  {
    const lf::base::AutoTimer t;
    for (long int i = 0; i < reps; i++) {
      auto res = getData_VEC(static_cast<double>(i));
      double s = 0.0;
      for (int j = 0; j < N; j++) {
        s += res[j];
      }
    }
  }

  std::cout << "II. Returning result through reference" << '\n';
  {
    const lf::base::AutoTimer t;
    std::vector<double> res(N);
    for (long int i = 0; i < reps; i++) {
      getData_REF(static_cast<double>(i), res);
      double s = 0.0;
      for (int j = 0; j < N; j++) {
        s += res[j];
      }
    }
  }

  std::cout << "III. Returning result through span" << '\n';
  {
    const lf::base::AutoTimer t;
    for (long int i = 0; i < reps; i++) {
      auto res = getData_SPAN(static_cast<double>(i));
      double s = 0.0;
      for (int j = 0; j < N; j++) {
        s += res[j];
      }
    }
  }

  return 0;
}
