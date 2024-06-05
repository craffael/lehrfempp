#include "print_info.h"

namespace lf::geometry {

void PrintInfo(std::ostream& o, const Geometry& geom, int output_ctrl) {
  const lf::base::RefEl geom_refel = geom.RefEl();
  o << "Reference element type for geometry: " << geom_refel << '\n';

  if (output_ctrl > 10) {
    // Dimensions and the derived type
    const base::dim_t dim_glob = geom.DimGlobal();
    const base::dim_t dim_local = geom.DimLocal();
    o << "world dim. = " << dim_glob << std::flush;
    o << ", loc dim =  " << dim_local << std::flush;
    o << "type: " << typeid(geom).name() << '\n';

    if (output_ctrl > 90) {
      // Geometry (coordinates of vertices)
      const Eigen::MatrixXd& ref_el_corners(geom_refel.NodeCoords());
      o << geom.Global(ref_el_corners) << '\n';
    }
  }

}  // void PrintInfo

std::ostream& operator<<(std::ostream& stream, const Geometry& geom) {
  return stream << geom.RefEl();
}

}  // namespace lf::geometry
