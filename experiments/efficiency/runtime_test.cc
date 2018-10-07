/** @file runtime_test.cc
 *  @brief Some tests for runtime behavior of certain C++ constructs
 */

#include <boost/timer/timer.hpp>
#include <iostream>
#include "lf/base/base.h"

static const int N = 10;

static double stat_tmp[N];

lf::base::RandomAccessRange<const double> getData_RAR(double offset) {
  for (int j = 0; j < N; j++) {
    stat_tmp[j] = 1.0 / (j + offset);
  }
  return {static_cast<double *>(stat_tmp), (stat_tmp + N)};
}

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

int main(int /*argc*/, const char * /*unused*/ []) {
  std::cout << "Runtime test for range access" << std::endl;

  std::cout << "I. Access through RandomAccessRange" << std::endl;
  const long int reps = 100000000L;
  {
    boost::timer::auto_cpu_timer t;
    for (long int i = 0; i < reps; i++) {
      auto res = getData_RAR(i);
      double s = 0.0;
      for (int j = 0; j < N; j++) {
        s += res[j];
      }
    }
  }

  std::cout << "II. Returning std::vector's" << std::endl;
  {
    boost::timer::auto_cpu_timer t;
    for (long int i = 0; i < reps; i++) {
      auto res = getData_VEC(i);
      double s = 0.0;
      for (int j = 0; j < N; j++) {
        s += res[j];
      }
    }
  }

  std::cout << "III. Returning result through reference" << std::endl;
  {
    boost::timer::auto_cpu_timer t;
    std::vector<double> res(N);
    for (long int i = 0; i < reps; i++) {
      getData_REF(i, res);
      double s = 0.0;
      for (int j = 0; j < N; j++) {
        s += res[j];
      }
    }
  }

  return 0;
}
