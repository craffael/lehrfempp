# include "outside.h"

using namespace lf::base;

void outside() {
  cc::Debug(1, "Outside: A level 1 debug message.");
  cc::Debug(5, "Outside: A level 5 debug message.");
  std::cout << "Outside: x = " << cv::Get("x") << "\n";
  std::cout << "Outside: y = " <<   cv::Get("y") << "\n";
  std::cout << "Outside: v1 = " <<   cv::Get("v1") << "\n";
}
