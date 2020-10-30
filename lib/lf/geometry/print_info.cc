#include "print_info.h"

namespace lf::geometry {

void PrintInfo(std::ostream& o, const Geometry& geom, int output_ctrl) {
  lf::base::RefEl geom_refel = geom.RefEl();
  o << "Reference element type for geometry: " << geom_refel << std::endl;

  if (output_ctrl > 10) {
    // Dimensions and the derived type
    int dim_glob = geom.DimGlobal();
    int dim_local = geom.DimLocal();
    o << "world dim. = " << dim_glob << std::flush;
    o << ", loc dim =  " << dim_local << std::flush;
    o << "type: " << typeid(geom).name() << std::endl;

    if (output_ctrl > 90) {
      // Geometry (coordinates of vertices)
      const Eigen::MatrixXd& ref_el_corners(geom_refel.NodeCoords());
      o << geom.Global(ref_el_corners) << std::endl;
    }
  }

}  // void PrintInfo

std::ostream& operator<<(std::ostream& stream, const Geometry& geom) {
  return stream << geom.RefEl();
}

}  // namespace lf::geometry
