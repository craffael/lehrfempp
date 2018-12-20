# include "outside.h"

using namespace lf::base;

void outside() {
  std::cout << "Outside: x = " << cv::Get<int>("x") << "\n";
  std::cout << "Outside: y = " <<   cv::Get<std::string>("y") << "\n";
}
