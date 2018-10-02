#include "quad_o1.h"
#include <cmath>
#include "point.h"
#include "segment_o1.h"
#include "tria_o1.h"

namespace lf::geometry {

bool assertNonDegenerateQuad(
    const Eigen::Matrix<double, Eigen::Dynamic, 4> &coords, double tol) {
  // World dimension
  const Geometry::dim_t wd = coords.rows();
  // Length tests
  double e0lensq = (coords.col(1) - coords.col(0)).squaredNorm();
  double e1lensq = (coords.col(2) - coords.col(1)).squaredNorm();
  double e2lensq = (coords.col(3) - coords.col(2)).squaredNorm();
  double e3lensq = (coords.col(0) - coords.col(3)).squaredNorm();
  // Test lengths of edges versus circumference.
  double circum = e0lensq + e1lensq + e2lensq + e3lensq;
  LF_VERIFY_MSG(e0lensq > tol * circum, "Collapsed edge 0");
  LF_VERIFY_MSG(e1lensq > tol * circum, "Collapsed edge 1");
  LF_VERIFY_MSG(e2lensq > tol * circum, "Collapsed edge 2");
  LF_VERIFY_MSG(e3lensq > tol * circum, "Collapsed edge 3");
  // Area test
  switch (wd) {
    case 2: {
      double ar1 =
          ((coords(0, 1) - coords(0, 0)) * (coords(1, 2) - coords(1, 0)) -
           (coords(1, 1) - coords(1, 0)) * (coords(0, 2) - coords(0, 0)));
      double ar2 =
          ((coords(0, 3) - coords(0, 0)) * (coords(1, 2) - coords(1, 0)) -
           (coords(1, 3) - coords(1, 0)) * (coords(0, 2) - coords(0, 0)));
      double area = std::fabs(ar1) + std::fabs(ar2);
      LF_VERIFY_MSG(area > tol * circum, "Degenerate 2D quad");
      return true;
      break;
    }
    case 3: {
      const Eigen::Matrix<double, 3, 4> c3d(coords.block<3, 4>(0, 0));
      double ar1 =
          ((c3d.col(1) - c3d.col(0)).cross(c3d.col(2) - c3d.col(0))).norm();
      double ar2 =
          ((c3d.col(3) - c3d.col(0)).cross(c3d.col(2) - c3d.col(0))).norm();
      double area = ar1 + ar2;
      LF_VERIFY_MSG(area > tol * circum, "Degenerate 3D quad");
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

QuadO1::QuadO1(Eigen::Matrix<double, Eigen::Dynamic, 4> coords)
    : coords_(std::move(coords)) {
  // Check validity of geometry (non-zero area)
  assertNonDegenerateQuad(coords_);
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
      result(i) = std::abs(jacobian.determinant());
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
    const RefinementPattern& ref_pat, base::dim_t codim) const {
  // The refinement pattern must be for a quadrilateral
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kQuad(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  // Allowed condimensions:
  // 0 -> child cells, 1-> child edges, 2 -> child points
  LF_VERIFY_MSG(codim < 3, "Illegal codim " << codim);

  // Lattice meshwidth
  const double h_lattice = 1.0 / static_cast<double>(ref_pat.LatticeConst());
  // Obtain geometry of children as lattice polygons
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
      child_polygons(ref_pat.ChildPolygons(codim));
  // Number of child segments
  const int no_children = child_polygons.size();
  // Check consistency of data
  LF_ASSERT_MSG(
      no_children == ref_pat.noChildren(codim),
      "no_children = " << no_children << " <-> " << ref_pat.noChildren(codim));

  // Variable for returning (unique pointers to) child geometries
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  // For each child entity create a geometry object and a unique pointer to it.
  for (int l = 0; l < no_children; l++) {
    // A single child entity is described by a lattice polygon with
    // a certain number of corners
    LF_VERIFY_MSG(
        child_polygons[l].rows() == 2,
        "child_polygons[" << l << "].rows() = " << child_polygons[l].rows());
    // Obtain physical (world) coordinates of the vertices of
    // of the child entities after normalizing lattice coordinates
    const Eigen::MatrixXd child_geo(
        Global(h_lattice * child_polygons[l].cast<double>()));
    // Treat child entities of different dimension
    switch (codim) {
      case 0: {  // Child cells
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
        break;
      }
      case 1: {
        // Child is an edge
        LF_VERIFY_MSG(
            child_polygons[l].cols() == 2,
            "child_polygons[l].cols() = " << child_polygons[l].cols());
        child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
        break;
      }
      case 2: {
        LF_VERIFY_MSG(
            child_polygons[l].cols() == 1,
            "child_polygons[l].cols() = " << child_polygons[l].cols());
        child_geo_uptrs.push_back(std::make_unique<Point>(child_geo));
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Illegal co-dimension");
        break;
      }
    }  // end switch codim
  }    // end loop over children
  return (child_geo_uptrs);
}  // end ChildGeometry()

//////////////////////////////////////////////////////////////////////
// Implementation class Parallelogram
//////////////////////////////////////////////////////////////////////

Parallelogram::Parallelogram(Eigen::Matrix<double, Eigen::Dynamic, 4> coords)
    : coords_(std::move(coords)),
      jacobian_(coords_.rows(), 2),
      jacobian_inverse_gramian_(coords_.rows(), 2),
      integrationElement_(0) {
  // Check validity of geometry (non-zero area)
  assertNonDegenerateQuad(coords_);
  // Check parallelogram property
  LF_ASSERT_MSG(
      (coords_.col(3) - coords_.col(2) + coords_.col(1) - coords_.col(0))
              .norm() < 1.0E-8 * (coords_.col(2) - coords_.col(0)).norm(),
      "No parallelogram!");
  init();
}

Parallelogram::Parallelogram(const Eigen::VectorXd& p0,
                             const Eigen::VectorXd& p1,
                             const Eigen::VectorXd& p2)
    : coords_(p0.rows(), 4),
      jacobian_(p0.rows(), 2),
      jacobian_inverse_gramian_(p0.rows(), 2),
      integrationElement_(0) {
  LF_ASSERT_MSG(p0.size() == p1.size(), "Vector length mismatch p0 <-> p1");
  LF_ASSERT_MSG(p0.size() == p2.size(), "Vector length mismatch p0 <-> p2");
  coords_.col(0) = p0;
  coords_.col(1) = p1;
  coords_.col(2) = p1 + (p2 - p0);
  coords_.col(3) = p2;
  assertNonDegenerateQuad(coords_);
  init();
}

void Parallelogram::init(void) {
  jacobian_ << coords_.col(1) - coords_.col(0), coords_.col(3) - coords_.col(0);
  // Distinguish between different world dimensions
  if (coords_.rows() == 2) {
    // 2D case: Simpler formula!
    jacobian_inverse_gramian_ = jacobian_.transpose().inverse();
    integrationElement_ = std::abs(jacobian_.determinant());
  } else {
    // 3D case: complicated formula
    jacobian_inverse_gramian_ = Eigen::MatrixXd(
        jacobian_ * (jacobian_.transpose() * jacobian_).inverse());
    integrationElement_ =
        std::sqrt((jacobian_.transpose() * jacobian_).determinant());
  }
}  // end init()

Eigen::MatrixXd Parallelogram::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(0) *
             (1 - local.array().row(0) - local.array().row(1)).matrix() +
         coords_.col(1) * local.row(0) + coords_.col(3) * local.row(1);
}

Eigen::MatrixXd Parallelogram::Jacobian(const Eigen::MatrixXd& local) const {
  return jacobian_.replicate(1, local.cols());
}

Eigen::MatrixXd Parallelogram::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {
  return jacobian_inverse_gramian_.replicate(1, local.cols());
}

Eigen::VectorXd Parallelogram::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  return Eigen::VectorXd::Constant(local.cols(), integrationElement_);
}


// essentially a copy of the same method for QuadO1
std::unique_ptr<Geometry> Parallelogram::SubGeometry(dim_t codim,
                                                     dim_t i) const {
  using std::make_unique;
  switch (codim) {
    case 0:
      LF_ASSERT_MSG(i == 0, "i is out of bounds.");
      return std::make_unique<Parallelogram>(coords_);
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
}  // end SubGeometry

// Essentially a copy of the same code for QuadO1
std::vector<std::unique_ptr<Geometry>> Parallelogram::ChildGeometry(
    const RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  // The refinement pattern must be for a quadrilateral
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kQuad(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  // Allowed condimensions:
  // 0 -> child cells, 1-> child edges, 2 -> child points
  LF_VERIFY_MSG(codim < 3, "Illegal codim " << codim);

  // Lattice meshwidth
  const double h_lattice = 1.0 / static_cast<double>(ref_pat.LatticeConst());
  // Obtain geometry of children as lattice polygons
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>
      child_polygons(ref_pat.ChildPolygons(codim));
  // Number of child segments
  const int no_children = child_polygons.size();
  // Check consistency of data
  LF_ASSERT_MSG(
      no_children == ref_pat.noChildren(codim),
      "no_children = " << no_children << " <-> " << ref_pat.noChildren(codim));

  // Variable for returning (unique pointers to) child geometries
  std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
  // For each child entity create a geometry object and a unique pointer to it.
  for (int l = 0; l < no_children; l++) {
    // A single child entity is described by a lattice polygon with
    // a certain number of corners
    LF_VERIFY_MSG(
        child_polygons[l].rows() == 2,
        "child_polygons[" << l << "].rows() = " << child_polygons[l].rows());
    // Obtain physical (world) coordinates of the vertices of
    // of the child entities after normalizing lattice coordinates
    const Eigen::MatrixXd child_geo(
        Global(h_lattice * child_polygons[l].cast<double>()));
    // Treat child entities of different dimension
    switch (codim) {
      case 0: {  // Child cells
        if (child_polygons[l].cols() == 3) {
          // Child cell is a triangle
          child_geo_uptrs.push_back(std::make_unique<TriaO1>(child_geo));
        } else if (child_polygons[l].cols() == 4) {
          // Child cell is a quadrilateral
          ////////////////////////////////////////////////////
          // Code specific for a parallleogram
          // If the lattice polygon is a parallelogram
          // then create another parallelogram, otherwise
          // a general quadrilateral
          if (isParallelogram(child_polygons[l])) {
            // Child cell is parallelogram as well
            child_geo_uptrs.push_back(
                std::make_unique<Parallelogram>(child_geo));
          } else {
            // General quadrilateral
            child_geo_uptrs.push_back(std::make_unique<QuadO1>(child_geo));
          }
        } else {
          LF_VERIFY_MSG(false, "child_polygons[" << l << "].cols() = "
                                                 << child_polygons[l].cols());
        }
        break;
      }
      case 1: {
        // Child is an edge
        LF_VERIFY_MSG(
            child_polygons[l].cols() == 2,
            "child_polygons[l].cols() = " << child_polygons[l].cols());
        child_geo_uptrs.push_back(std::make_unique<SegmentO1>(child_geo));
        break;
      }
      case 2: {
        LF_VERIFY_MSG(
            child_polygons[l].cols() == 1,
            "child_polygons[l].cols() = " << child_polygons[l].cols());
        child_geo_uptrs.push_back(std::make_unique<Point>(child_geo));
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Illegal co-dimension");
        break;
      }
    }  // end switch codim
  }    // end loop over children
  return (child_geo_uptrs);
}  // end ChildGeometry()

}  // namespace lf::geometry
