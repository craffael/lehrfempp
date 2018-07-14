#include "segment_o1.h"
#include "point.h"

namespace lf::geometry {

Eigen::MatrixXd SegmentO1::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(1) * local + coords_.col(0) * (1 - local.array()).matrix();
}

Eigen::MatrixXd SegmentO1::Jacobian(const Eigen::MatrixXd& local) const {
  return (coords_.col(1) - coords_.col(0)).replicate(1, local.cols());
}

Eigen::MatrixXd SegmentO1::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {
  if (DimGlobal() == 1) {
    return (coords_.col(1) - coords_.col(0))
        .cwiseInverse()
        .replicate(1, local.cols());
  }
  return ((coords_.col(1) - coords_.col(0)) /
          (coords_.col(1) - coords_.col(0)).squaredNorm())
      .replicate(1, local.cols());
}

Eigen::VectorXd SegmentO1::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  return Eigen::VectorXd::Constant(local.cols(),
                                   (coords_.col(1) - coords_.col(0)).norm());
}

std::unique_ptr<Geometry> SegmentO1::SubGeometry(dim_t codim, dim_t i) const {
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "i is out of bounds.");
    return std::make_unique<SegmentO1>(coords_);
  }
  if (codim == 1) {
    LF_ASSERT_MSG(i >= 0 && i < 2, "i is out of bounds.");
    return std::make_unique<Point>(coords_.col(i));
  }
  LF_VERIFY_MSG(false, "codim is out of bounds.");
}

std::vector<std::unique_ptr<Geometry>> SegmentO1::ChildGeometry(
    const RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  // The refinement pattern must be for a segment
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kSegment(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  // Lattice meshwidth
  const double h_lattice = 1.0 / (double)ref_pat.LatticeConst();
  // Obtain geometry of children as vector of pairs of lattice coordinates
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
      child_polygons(ref_pat.ChildPolygons(0));
  // Number of child segments
  const int no_children = child_polygons.size();
  LF_VERIFY_MSG(
      no_children == ref_pat.noChildren(0),
      "no_children = " << no_children << " <-> " << ref_pat.noChildren(0));
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  // For each child create a geometry object and a unique pointer to it.
  for (int l = 0; l < no_children; l++) {
    // A single child must be described by an interval, that is
    // two differente lattice coordinates
    LF_VERIFY_MSG(child_polygons[l].rows() == 1,
                  "child_polygons[l].rows() = " << child_polygons[l].rows());
    LF_VERIFY_MSG(child_polygons[l].cols() == 2,
                  "child_polygons[l].cols() = " << child_polygons[l].cols());
    // Normalize lattice coordinates
    const Eigen::MatrixXd child_geo(
        Global(h_lattice * child_polygons[l].cast<double>()));
    child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
  }
  return child_geo_uptrs;
}

/* OLD IMPLEMENTATION based on explicit refinement patterns
std::vector<std::unique_ptr<Geometry>>
SegmentO1::ChildGeometry(const RefinementPattern &ref_pat) const {
  RefPat ref_pattern = ref_pat.refpat();
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  switch (ref_pattern) {
  case (int)RefPat::rp_copy: {
    child_geo_uptrs.push_back(std::make_unique<SegmentO1>(coords_));
    break;
  }
  case (int)RefPat::rp_regular:
  case (int)RefPat::rp_split: {
    const dim_t dim_global(coords_.rows());
    Eigen::Matrix<double,Eigen::Dynamic,2> child_geo(dim_global,2);
    // Reference coordinates of endpoints and midpoint
    Eigen::Matrix<double,1,3> ref_coords; ref_coords << 0.0,0.5,1.0;
    // Fetch point coordinates
    Eigen::MatrixXd point_coords(Global(ref_coords));
    // "first half edge" from endpoint 0 to midpoint
    child_geo = point_coords.block(0,0,dim_global,2);
    child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
    // "second half edge" from midpoint to endpoint 1
    child_geo = point_coords.block(0,1,dim_global,2);
    child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
    break;
  }
  default: {
    LF_VERIFY_MSG(false,"Invalid refinement pattern " << (int)ref_pattern);
    break;
  }
  } // end switch
  return std::move(child_geo_uptrs);
  } */

}  // namespace lf::geometry
