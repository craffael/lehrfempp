#include "tria_o1.h"

#include <Eigen/Eigen>

#include "point.h"
#include "segment_o1.h"

namespace lf::geometry {

bool assertNonDegenerateTriangle(
    const Eigen::Matrix<double, Eigen::Dynamic, 3>& coords, double tol) {
  // World dimension
  const Geometry::dim_t wd = coords.rows();
  // Length tests
  double e0lensq = (coords.col(1) - coords.col(0)).squaredNorm();
  double e1lensq = (coords.col(2) - coords.col(1)).squaredNorm();
  double e2lensq = (coords.col(0) - coords.col(2)).squaredNorm();
  // Test lengths of edges versus circumference.
  double circum = e0lensq + e1lensq + e2lensq;
  LF_VERIFY_MSG(e0lensq > tol * circum, "Collapsed edge 0");
  LF_VERIFY_MSG(e1lensq > tol * circum, "Collapsed edge 1");
  LF_VERIFY_MSG(e2lensq > tol * circum, "Collapsed edge 2");
  // Area test
  switch (wd) {
    case 2: {
      double area = std::fabs(
          ((coords(0, 1) - coords(0, 0)) * (coords(1, 2) - coords(1, 0)) -
           (coords(1, 1) - coords(1, 0)) * (coords(0, 2) - coords(0, 0))));
      LF_VERIFY_MSG(area > tol * circum,
                    "Degenerate 2D triangle, geo = " << coords);
      return true;
      break;
    }
    case 3: {
      const Eigen::Matrix<double, 3, 3> c3d(coords.block<3, 3>(0, 0));
      double area =
          ((c3d.col(1) - c3d.col(0)).cross(c3d.col(2) - c3d.col(0))).norm();
      LF_VERIFY_MSG(area > tol * circum, "Degenerate 3D triangle");
      return true;
      break;
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal world dimension" << wd);
      break;
    }
  }
  return false;
}

TriaO1::TriaO1(Eigen::Matrix<double, Eigen::Dynamic, 3> coords)
    : coords_(std::move(coords)),
      jacobian_(coords_.rows(), 2),
      jacobian_inverse_gramian_(coords_.rows(), 2) {
  // Make sure that the triangle has a proper shape.
  assertNonDegenerateTriangle(coords_);
  // Precompute constant Jacobian
  jacobian_ << coords_.col(1) - coords_.col(0), coords_.col(2) - coords_.col(0);
  // Precompute constant metric factor
  if (coords_.rows() == 2) {
    jacobian_inverse_gramian_ = jacobian_.transpose().inverse();
    integrationElement_ = std::abs(jacobian_.determinant());
  } else {
    jacobian_inverse_gramian_ = Eigen::MatrixXd(
        jacobian_ * (jacobian_.transpose() * jacobian_).inverse());
    integrationElement_ =
        std::sqrt((jacobian_.transpose() * jacobian_).determinant());
  }
}  // end constructor

Eigen::MatrixXd TriaO1::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(0) *
             (1 - local.array().row(0) - local.array().row(1)).matrix() +
         coords_.col(1) * local.row(0) + coords_.col(2) * local.row(1);
}

std::unique_ptr<Geometry> TriaO1::SubGeometry(dim_t codim, dim_t i) const {
  using std::make_unique;
  switch (codim) {
    case 0:
      LF_ASSERT_MSG(i == 0, "codim = 0: i = " << i << " is out of bounds.");
      return std::make_unique<TriaO1>(coords_);
    case 1:
      LF_ASSERT_MSG(i >= 0 && i < 3,
                    "codim =1: i = " << i << " is out of bounds.");
      return make_unique<SegmentO1>(
          (Eigen::Matrix<double, Eigen::Dynamic, 2>(DimGlobal(), 2)
               << coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 0)),
           coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 1)))
              .finished());
    case 2:
      LF_ASSERT_MSG(i >= 0 && i < 3,
                    "codim = 2: i = " << i << " is out of bounds.");
      return make_unique<Point>(coords_.col(i));
    default:
      LF_VERIFY_MSG(false, "codim " << codim << " is out of bounds.");
  }
}

std::vector<std::unique_ptr<Geometry>> TriaO1::ChildGeometry(
    const RefinementPattern& ref_pat, base::dim_t codim) const {
  // The refinement pattern must be for a triangle
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kTria(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  LF_VERIFY_MSG(codim < 3, "Illegal codim " << codim);
  // Lattice meshwidth
  const double h_lattice = 1.0 / static_cast<double>(ref_pat.LatticeConst());
  // Obtain geometry of children as lattice polygon
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
      child_polygons(ref_pat.ChildPolygons(codim));
  // Number of child segments
  const auto no_children = base::narrow<base::size_type>(child_polygons.size());
  LF_VERIFY_MSG(
      no_children == ref_pat.NumChildren(codim),
      "no_children = " << no_children << " <-> " << ref_pat.NumChildren(codim));
  // return variable
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  // For each child triangle create a geometry object and a unique pointer to
  // it.
  for (int l = 0; l < no_children; l++) {
    // A single child triangle is described by a lattice polygon with
    // three vertices, a child segment by a polygon with two vertices,
    // a child point by a single point = "polygon with one corner"
    LF_VERIFY_MSG(child_polygons[l].rows() == 2,
                  "child_polygons[l].rows() = " << child_polygons[l].rows());
    LF_VERIFY_MSG(child_polygons[l].cols() == 3 - codim,
                  "child_polygons[l].cols() = " << child_polygons[l].cols());
    // Normalize lattice coordinates
    const Eigen::MatrixXd child_geo(
        Global(h_lattice * child_polygons[l].cast<double>()));
    switch (codim) {
      case 0: {
        // child is a triangle
        child_geo_uptrs.push_back(std::make_unique<TriaO1>(child_geo));
        break;
      }
      case 1: {
        // child is an edge
        child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
        break;
      }
      case 2: {
        // child is a node
        child_geo_uptrs.push_back(std::make_unique<Point>(child_geo));
        break;
      }
      default:
        LF_VERIFY_MSG(false, "codim out of bounds.");
    }  // end switch codim
  }    // end loop over the children
  return (child_geo_uptrs);
}

}  // namespace lf::geometry
