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
    const Eigen::MatrixXd& local) const {
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
    const RefinementPattern& ref_pat, base::dim_t codim) const {
  // The refinement pattern must be for a segment
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kSegment(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  LF_VERIFY_MSG((codim < 2), "Illegal codim = " << codim);

  // Lattice meshwidth
  const double h_lattice = 1.0 / static_cast<double>(ref_pat.LatticeConst());
  // Return variable
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};

  // Obtain geometry of children as vector of pairs of lattice coordinates
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
      child_polygons(ref_pat.ChildPolygons(codim));
  // Number of child entities
  const int no_children = child_polygons.size();
  LF_VERIFY_MSG(
      no_children == ref_pat.noChildren(codim),
      "no_children = " << no_children << " <-> " << ref_pat.noChildren(codim));
  // For each child create a geometry object and a unique pointer to it.
  for (int l = 0; l < no_children; l++) {
    // codim == 0:A single child must be described by an interval, that is
    // two different lattice coordinates
    // codim == 1: a point's location is just one number
    LF_VERIFY_MSG(child_polygons[l].rows() == 1,
                  "child_polygons[l].rows() = " << child_polygons[l].rows());
    LF_VERIFY_MSG(child_polygons[l].cols() == (2 - codim),
                  "child_polygons[l].cols() = " << child_polygons[l].cols());
    // Normalize lattice coordinates
    const Eigen::MatrixXd child_geo(
        Global(h_lattice * child_polygons[l].cast<double>()));

    switch (codim) {
      case 0: {
        child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
        break;
      }
      case 1: {
        child_geo_uptrs.push_back(std::make_unique<Point>(child_geo));
        break;
      }
      default:
        LF_VERIFY_MSG(false, "codim out of bounds.");
    }  // end switch codim
  }    // end loop over child entities
  return child_geo_uptrs;
}

}  // namespace lf::geometry
