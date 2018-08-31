/** @file runtime_test.cc
 *  @brief Some tests for runtime behavior of certain C++ constructs
 */

#include <iostream>
#include "lf/base/base.h"
#include <boost/timer/timer.hpp>

static const int N = 10;

static double stat_tmp[N];

lf::base::RandomAccessRange<const double> getData_RAR(double offset) {
  for (int j=0;j<N; j++) {
    stat_tmp[j] = 1.0/(j+offset);
  }
  return {(const double *)stat_tmp,(const double *)(stat_tmp+N)};
}

std::vector<double> getData_VEC(double offset) {
  std::vector<double> tmp(10);
  for (int j=0;j<N; j++) {
    tmp[j] = 1.0/(j+offset);
  }
  return tmp;
}

void getData_REF(double offset,std::vector<double> &res) {
  for (int j=0;j<N; j++) {
    res[j] = 1.0/(j+offset);
  }
}

int main(int /*argc*/, const char* /*unused*/ []) {
  std::cout << "Runtime test" << std::endl;

  {
    boost::timer::auto_cpu_timer t;
    for (long i=0;i<10000;i++) {
      auto res = getData_RAR((double)i);
      double s=0.0; for (int j=0;j<N;j++) { s += res[j]; }
    }
  }
  
  {
    boost::timer::auto_cpu_timer t;
    for (long i=0;i<10000;i++) {
      auto res = getData_VEC((double)i);
      double s=0.0; for (int j=0;j<N;j++) { s += res[j]; }
    }
  }

  {
    boost::timer::auto_cpu_timer t;
    std::vector<double> res(N);
    for (long i=0;i<10000;i++) {
      getData_REF((double)i,res);
      double s=0.0; for (int j=0;j<N;j++) { s += res[j]; }
    }
  }
  
  return 0;
}
