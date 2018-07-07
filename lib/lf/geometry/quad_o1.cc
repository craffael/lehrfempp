#include "quad_o1.h"
#include "point.h"
#include "segment_o1.h"
#include <cmath>

namespace lf::geometry {

  QuadO1::QuadO1(Eigen::Matrix<double, Eigen::Dynamic, 4> coords):
    coords_(std::move(coords)) {
    // Check validity of geometry (non-zero area)
    double ar1 = ((coords_(0,1)-coords_(0,0))*(coords_(1,2)-coords_(1,0))
		  - (coords_(1,1)-coords_(1,0))*(coords_(0,2)-coords_(0,0)));
    double ar2 = ((coords_(0,3)-coords_(0,0))*(coords_(1,2)-coords_(1,0))
		  - (coords_(1,3)-coords_(1,0))*(coords_(0,2)-coords_(0,0)));
    double area = std::fabs(ar1)+std::fabs(ar2);
    double e0lensq = (std::pow(coords_(0,1)-coords_(0,0),2)+
		      std::pow(coords_(1,1)-coords_(1,0),2));
    double e1lensq = (std::pow(coords_(0,2)-coords_(0,1),2)+
		      std::pow(coords_(1,2)-coords_(1,1),2));
    double e2lensq = (std::pow(coords_(0,3)-coords_(0,2),2)+
		      std::pow(coords_(1,3)-coords_(1,2),2));
    double e3lensq = (std::pow(coords_(0,0)-coords_(0,3),2)+
		      std::pow(coords_(1,0)-coords_(1,3),2));
    double circum = e0lensq + e1lensq + e2lensq + e3lensq;
    LF_VERIFY_MSG(e0lensq > 1.0E-8 * circum,"Collapsed edge 0");
    LF_VERIFY_MSG(e1lensq > 1.0E-8 * circum,"Collapsed edge 1");
    LF_VERIFY_MSG(e2lensq > 1.0E-8 * circum,"Collapsed edge 2");
    LF_VERIFY_MSG(e3lensq > 1.0E-8 * circum,"Collapsed edge 3");
    LF_VERIFY_MSG(area > 1.0E-8*circum,"Degenerate quad");
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

  std::vector<std::unique_ptr<Geometry>>
  QuadO1::ChildGeometry(const RefinementPattern &ref_pat) const {
    LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kQuad(),
		  "Refinement pattern not for triangle")
    RefPat ref_pattern = ref_pat.refpat();
    const int anchor = ref_pat.anchor();
    // Vector for returning unique pointers to child geometries
    std::vector<std::unique_ptr<Geometry>> child_geo_uptrs{};
    // Key properties of the current cell
    const dim_t dim_global(DimGlobal());
    const lf::base::RefEl ref_el = RefEl();
    // Coordinates of corners and midpoints in reference quadrilateral
    Eigen::Matrix<double,2,4> ref_corner_coords(ref_el.NodeCoords());
    Eigen::Matrix<double,2,4> ref_midpoint_coords(2,4);
    ref_midpoint_coords << 0.5,1.0,0.5,0.0,0.0,0.5,1.0,0.5;

    // Physical coordinates of corners and midpoints
    Eigen::MatrixXd corner_coords = Global(ref_corner_coords);
    Eigen::MatrixXd midpoint_coords = Global(ref_midpoint_coords);
    Eigen::MatrixXd tria_child_coords(2,3);
    Eigen::MatrixXd quad_child_coords(2,4);

    // Remap local indices according to anchor values
    const int mod_0 = (0+anchor)%4;
    const int mod_1 = (1+anchor)%4;
    const int mod_2 = (2+anchor)%4;
    const int mod_3 = (3+anchor)%4;
    
    // Create child geometries according to refinement patterns and selection
    switch (ref_pattern) {
    case (int)RefPat::rp_nil: {
      break;
    }
    case (int)RefPat::rp_copy: {
      child_geo_uptrs.push_back(std::make_unique<QuadO1>(coords_));
      break;
    }
    case (int)RefPat::rp_trisect: {
      // Partition a quad into three triangle, the anchor edge
      // being split in the process
	tria_child_coords.col(0) = midpoint_coords.col(mod_0);
	tria_child_coords.col(1) = corner_coords.col(mod_2);
	tria_child_coords.col(2) = corner_coords.col(mod_3);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = midpoint_coords.col(mod_0);
	tria_child_coords.col(1) = corner_coords.col(mod_0);
	tria_child_coords.col(2) = corner_coords.col(mod_3);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = midpoint_coords.col(mod_0);
	tria_child_coords.col(1) = corner_coords.col(mod_1);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));
	tria_child_coords.col(2) = corner_coords.col(mod_2);
      break;
    }
    case (int)RefPat::rp_quadsect: {
      // Partition a quad into four triangle, thus 
      // splitting two edges. The one with the smaller sub index is the
      // anchor edge
	tria_child_coords.col(0) = corner_coords.col(mod_0);
	tria_child_coords.col(1) = corner_coords.col(mod_3);
	tria_child_coords.col(2) = midpoint_coords.col(mod_0);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = corner_coords.col(mod_1);
	tria_child_coords.col(1) = midpoint_coords.col(mod_1);
	tria_child_coords.col(2) = midpoint_coords.col(mod_0);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = corner_coords.col(mod_2);
	tria_child_coords.col(1) = corner_coords.col(mod_3);
	tria_child_coords.col(2) = midpoint_coords.col(mod_1);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = midpoint_coords.col(mod_0);
	tria_child_coords.col(1) = midpoint_coords.col(mod_2);
	tria_child_coords.col(2) = corner_coords.col(mod_3);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));
      break;
    }
    case (int)RefPat::rp_bisect: 
    case (int)RefPat::rp_split: {
      // Cut a quadrilateral into two 
	quad_child_coords.col(0) = corner_coords.col(mod_0);
	quad_child_coords.col(1) = midpoint_coords.col(mod_0);
	quad_child_coords.col(2) = midpoint_coords.col(mod_2);
	quad_child_coords.col(3) = corner_coords.col(mod_3);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));

	quad_child_coords.col(0) = corner_coords.col(mod_1);
	quad_child_coords.col(1) = corner_coords.col(mod_2);
	quad_child_coords.col(2) = midpoint_coords.col(mod_2);
	quad_child_coords.col(3) = midpoint_coords.col(mod_0);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));
      break;
    }
    case (int)RefPat::rp_threeedge: {
	quad_child_coords.col(0) = corner_coords.col(mod_2);
	quad_child_coords.col(1) = corner_coords.col(mod_3);
	quad_child_coords.col(2) = midpoint_coords.col(mod_3);
	quad_child_coords.col(3) = midpoint_coords.col(mod_1);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));

	tria_child_coords.col(0) = corner_coords.col(mod_0);
	tria_child_coords.col(1) = midpoint_coords.col(mod_0);
	tria_child_coords.col(2) = midpoint_coords.col(mod_3);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = corner_coords.col(mod_1);
	tria_child_coords.col(1) = midpoint_coords.col(mod_0);
	tria_child_coords.col(2) = midpoint_coords.col(mod_1);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));

	tria_child_coords.col(0) = midpoint_coords.col(mod_0);
	tria_child_coords.col(1) = midpoint_coords.col(mod_1);
	tria_child_coords.col(2) = midpoint_coords.col(mod_3);
	child_geo_uptrs.push_back(std::make_unique<TriaO1>(tria_child_coords));
      break;
    }
    case (int)RefPat::rp_barycentric:
    case (int)RefPat::rp_regular: {
      // Fully symmetric splitting into four quadrilaterals
      // Obtain coordinates of center of gravity
      Eigen::Matrix<double,2,1> ref_baryc_coords = Eigen::Vector2d({0.5,0.5});
      Eigen::MatrixXd baryc_coords = Global(ref_baryc_coords);

      quad_child_coords.col(0) = corner_coords.col(0);
	quad_child_coords.col(1) = midpoint_coords.col(0);
	quad_child_coords.col(2) = baryc_coords;
	quad_child_coords.col(3) = midpoint_coords.col(3);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));

	quad_child_coords.col(0) = corner_coords.col(1);
	quad_child_coords.col(1) = midpoint_coords.col(1);
	quad_child_coords.col(2) = baryc_coords;
	quad_child_coords.col(3) = midpoint_coords.col(0);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));

	quad_child_coords.col(0) = corner_coords.col(2);
	quad_child_coords.col(1) = midpoint_coords.col(1);
	quad_child_coords.col(2) = baryc_coords;
	quad_child_coords.col(3) = midpoint_coords.col(2);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));

	quad_child_coords.col(0) = corner_coords.col(3);
	quad_child_coords.col(1) = midpoint_coords.col(2);
	quad_child_coords.col(2) = baryc_coords;
	quad_child_coords.col(3) = midpoint_coords.col(3);
	child_geo_uptrs.push_back(std::make_unique<QuadO1>(quad_child_coords));
      break;
    }
    default: {
      LF_VERIFY_MSG(false,"No valid refinement pattern for a quad");
      break;
    }
    } // end switch ref_pattern
    return std::move(child_geo_uptrs);
  }
    
}  // namespace lf::geometry
