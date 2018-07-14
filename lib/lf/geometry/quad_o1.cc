#include "quad_o1.h"
#include <cmath>
#include "point.h"
#include "segment_o1.h"
#include "tria_o1.h"

namespace lf::geometry {

QuadO1::QuadO1(Eigen::Matrix<double, Eigen::Dynamic, 4> coords)
    : coords_(std::move(coords)) {
  // Check validity of geometry (non-zero area)
  double ar1 =
      ((coords_(0, 1) - coords_(0, 0)) * (coords_(1, 2) - coords_(1, 0)) -
       (coords_(1, 1) - coords_(1, 0)) * (coords_(0, 2) - coords_(0, 0)));
  double ar2 =
      ((coords_(0, 3) - coords_(0, 0)) * (coords_(1, 2) - coords_(1, 0)) -
       (coords_(1, 3) - coords_(1, 0)) * (coords_(0, 2) - coords_(0, 0)));
  double area = std::fabs(ar1) + std::fabs(ar2);
  double e0lensq = (std::pow(coords_(0, 1) - coords_(0, 0), 2) +
                    std::pow(coords_(1, 1) - coords_(1, 0), 2));
  double e1lensq = (std::pow(coords_(0, 2) - coords_(0, 1), 2) +
                    std::pow(coords_(1, 2) - coords_(1, 1), 2));
  double e2lensq = (std::pow(coords_(0, 3) - coords_(0, 2), 2) +
                    std::pow(coords_(1, 3) - coords_(1, 2), 2));
  double e3lensq = (std::pow(coords_(0, 0) - coords_(0, 3), 2) +
                    std::pow(coords_(1, 0) - coords_(1, 3), 2));
  double circum = e0lensq + e1lensq + e2lensq + e3lensq;
  LF_VERIFY_MSG(e0lensq > 1.0E-8 * circum, "Collapsed edge 0");
  LF_VERIFY_MSG(e1lensq > 1.0E-8 * circum, "Collapsed edge 1");
  LF_VERIFY_MSG(e2lensq > 1.0E-8 * circum, "Collapsed edge 2");
  LF_VERIFY_MSG(e3lensq > 1.0E-8 * circum, "Collapsed edge 3");
  LF_VERIFY_MSG(area > 1.0E-8 * circum, "Degenerate quad");
}

Eigen::MatrixXd QuadO1::Global(const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 2, "reference coords must be 2-vectors");
  return coords_.col(0) *
             ((1 - local.array().row(0)) * (1 - local.array().row(1)))
                 .matrix() +
         coords_.col(1) *
             (local.array().row(0) * (1 - local.array().row(1))).matrix() +
         coords_.col(2) *
             (local.array().row(0) * local.array().row(1)).matrix() +
         coords_.col(3) *
             ((1 - local.array().row(0)) * local.array().row(1)).matrix();
}

Eigen::MatrixXd QuadO1::Jacobian(const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(DimGlobal(), local.cols() * 2);

  for (int i = 0; i < local.cols(); ++i) {
    result.col(2 * i) = (coords_.col(1) - coords_.col(0)) * (1 - local(1, i)) +
                        (coords_.col(2) - coords_.col(3)) * local(1, i);
    result.col(2 * i + 1) =
        (coords_.col(3) - coords_.col(0)) * (1 - local(0, i)) +
        (coords_.col(2) - coords_.col(1)) * local(0, i);
  }
  return result;
}

Eigen::MatrixXd QuadO1::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(DimGlobal(), local.cols() * 2);
  Eigen::MatrixXd jacobian(DimGlobal(), 2);

  for (int i = 0; i < local.cols(); ++i) {
    jacobian.col(0) = (coords_.col(1) - coords_.col(0)) * (1 - local(1, i)) +
                      (coords_.col(2) - coords_.col(3)) * local(1, i);
    jacobian.col(1) = (coords_.col(3) - coords_.col(0)) * (1 - local(0, i)) +
                      (coords_.col(2) - coords_.col(1)) * local(0, i);

    if (DimGlobal() == 2) {
      result.block(0, 2 * i, DimGlobal(), 2) = jacobian.transpose().inverse();
    } else {
      result.block(0, 2 * i, DimGlobal(), 2) = (jacobian.transpose() * jacobian)
                                                   .colPivHouseholderQr()
                                                   .solve(jacobian.transpose())
                                                   .transpose();
    }
  }

  return result;
}

Eigen::VectorXd QuadO1::IntegrationElement(const Eigen::MatrixXd& local) const {
  Eigen::VectorXd result(local.cols());
  Eigen::MatrixXd jacobian(DimGlobal(), 2);
  for (int i = 0; i < local.cols(); ++i) {
    jacobian.col(0) = (coords_.col(1) - coords_.col(0)) * (1 - local(1, i)) +
                      (coords_.col(2) - coords_.col(3)) * local(1, i);
    jacobian.col(1) = (coords_.col(3) - coords_.col(0)) * (1 - local(0, i)) +
                      (coords_.col(2) - coords_.col(1)) * local(0, i);

    if (DimGlobal() == 2) {
      result(i) = jacobian.determinant();
    } else {
      result(i) = std::sqrt((jacobian.transpose() * jacobian).determinant());
    }
  }
  return result;
}

std::unique_ptr<Geometry> QuadO1::SubGeometry(dim_t codim, dim_t i) const {
  using std::make_unique;
  switch (codim) {
    case 0:
      LF_ASSERT_MSG(i == 0, "i is out of bounds.");
      return std::make_unique<QuadO1>(coords_);
    case 1:
      LF_ASSERT_MSG(i >= 0 && i < 4, "i is out of bounds.");
      return make_unique<SegmentO1>(
          (Eigen::Matrix<double, Eigen::Dynamic, 2>(DimGlobal(), 2)
               << coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 0)),
           coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 1)))
              .finished());
    case 2:
      LF_ASSERT_MSG(i >= 0 && i < 4, "i is out of bounds.");
      return make_unique<Point>(coords_.col(i));
    default:
      LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}

std::vector<std::unique_ptr<Geometry>> QuadO1::ChildGeometry(
    const RefinementPattern& ref_pat) const {
  // The refinement pattern must be for a quadrilateral
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kQuad(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  // Lattice meshwidth
  const double h_lattice = 1.0 / static_cast<double>(ref_pat.LatticeConst());
  // Obtain geometry of children as lattice polygons
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
      child_polygons(ref_pat.ChildPolygons(0));
  // Number of child segments
  const int no_children = child_polygons.size();
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  // For each child cell create a geometry object and a unique pointer to it.
  for (int l = 0; l < no_children; l++) {
    // A single child cell is described by a lattice polygon with
    // three or four vertices
    LF_VERIFY_MSG(
        child_polygons[l].rows() == 2,
        "child_polygons[" << l << "].rows() = " << child_polygons[l].rows());

    // Normalize lattice coordinates
    const Eigen::MatrixXd child_geo(
        Global(h_lattice * child_polygons[l].cast<double>()));
    if (child_polygons[l].cols() == 3) {
      // Child cell is a triangle
      child_geo_uptrs.push_back(std::make_unique<TriaO1>(child_geo));
    } else if (child_polygons[l].cols() == 4) {
      // Child cell is a quadrilateral
      child_geo_uptrs.push_back(std::make_unique<QuadO1>(child_geo));
    } else {
      LF_VERIFY_MSG(false, "child_polygons[" << l << "].cols() = "
                                             << child_polygons[l].cols());
    }
  }
  return (child_geo_uptrs);
}  // end ChildGeometry()

}  // namespace lf::geometry
