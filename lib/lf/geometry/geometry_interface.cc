/* @file geometry_interface.cc
 *  @brief implementation of functions operating on generic geometry objects
 */

#include "geometry_interface.h"

namespace lf::geometry {

std::ostream& operator<<(std::ostream& stream, const Geometry& geom) {
  if (Geometry::output_ctrl_ == 0) {
    return stream << geom.RefEl();
  } else {
    PrintInfo(geom, stream);
  }
}

double Volume(const Geometry& geo) {
  const lf::base::dim_t dim = geo.DimGlobal();
  const lf::base::dim_t refdim = geo.DimLocal();

  Eigen::Matrix<double, Eigen::Dynamic, 1> refc(refdim, 1);
  double refvol = 1.0;
  switch (geo.RefEl()) {
    case lf::base::RefEl::kPoint(): {
      return 1.0;
    }
    case lf::base::RefEl::kSegment(): {
      refc << 0.5;
      break;
    }
    case lf::base::RefEl::kTria(): {
      refvol = 0.5;
      refc << 1.0 / 3, 1.0 / 3;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      refc << 0.5, 0.5;
      break;
    }
    default: {
      LF_VERIFY_MSG(false,
                    "Volume not available for " << geo.RefEl().ToString());
    }
  }  // end switch RefEl
  return (refvol * ((geo.IntegrationElement(refc))[0]));
}

}  // namespace lf::geometry
