#include "print_info.h"

namespace lf::geometry {

CONTROLDECLARECOMMENT(Geometry, output_ctrl_, "output_ctrl_",
                      "Diagnostics control for Geometry");

void PrintInfo(const Geometry &geom, std::ostream &o) {
  int dim_glob = geom.DimGlobal();
  int dim_local = geom.DimLocal();

  o << "Global dimension: " << dim_glob << std::endl;
  o << "Local dimension: " << dim_local << std::endl;
  o << "Type of reference element: " << geom.RefEl() << std::endl;
  // Working
}

}  // namespace lf::geometry
