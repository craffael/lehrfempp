#include "tria_o1.h"
#include "point.h"
#include "segment_o1.h"

namespace lf::geometry {

TriaO1::TriaO1(Eigen::Matrix<double, Eigen::Dynamic, 3> coords)
    : coords_(std::move(coords)),
      jacobian_(coords_.rows(), 2),
      jacobian_inverse_gramian_(coords_.rows(), 2),
      integrationElement_(0) {
  jacobian_ << coords_.col(1) - coords_.col(0), coords_.col(2) - coords_.col(0);

  if (coords_.rows() == 2) {
    jacobian_inverse_gramian_ = jacobian_.transpose().inverse();
    integrationElement_ = jacobian_.determinant();
  } else {
    jacobian_inverse_gramian_ = Eigen::MatrixXd(
        jacobian_ * (jacobian_.transpose() * jacobian_).inverse());
    integrationElement_ =
        std::sqrt((jacobian_.transpose() * jacobian_).determinant());
  }
}

Eigen::MatrixXd TriaO1::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(0) *
             (1 - local.array().row(0) - local.array().row(1)).matrix() +
         coords_.col(1) * local.row(0) + coords_.col(2) * local.row(1);
}

std::unique_ptr<Geometry> TriaO1::SubGeometry(dim_t codim, dim_t i) const {
  using std::make_unique;
  switch (codim) {
    case 0:
      LF_ASSERT_MSG(i == 0, "i is out of bounds.");
      return std::make_unique<TriaO1>(coords_);
    case 1:
      LF_ASSERT_MSG(i >= 0 && i < 3, "i is out of bounds.");
      return make_unique<SegmentO1>(
          (Eigen::Matrix<double, Eigen::Dynamic, 2>(DimGlobal(), 2)
               << coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 0)),
           coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 1)))
              .finished());
    case 2:
      LF_ASSERT_MSG(i >= 0 && i < 3, "i is out of bounds.");
      return make_unique<Point>(coords_.col(i));
    default:
      LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}

  std::unique_ptr<Geometry>
  TriaO1::ChildGeometry(int ref_pattern,int selector) const {
    const dim_t dim_global(DimGlobal());
    const lf::base::RefEl ref_el = RefEl();
    // Coordinates of corners and midpoints in reference triangle
    Eigen::MatrixXd ref_corner_coords(ref_el.NodeCoords());
    Eigen::MatrixXd ref_midpoint_coords(2,3);
    ref_midpoint_coords << 0.5,0.5,0.0,0.0,0.5,0.5;

    // Physical coordinates of corners and midpoints
    Eigen::MatrixXd corner_coords = Global(ref_corner_coords);
    Eigen::MatrixXd midpoint_coords = Global(ref_midpoint_coords);
    Eigen::MatrixXd child_coords(2,3);

    // Create child geometries according to refinement patterns and selection
    switch (ref_pattern) {
    case (int)RefinementPattern::rp_copy: {
      return std::make_unique<TriaO1>(coords_);
      break;
    }
    case (int)RefinementPattern::rp_bisect_0: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_bisect_1: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = corner_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = corner_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_bisect_2: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = corner_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = corner_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_trisect_01: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_trisect_02: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_trisect_10: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = corner_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_trisect_12: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = corner_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_trisect_20: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_trisect_21: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = corner_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_quadsect_0: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 3: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_quadsect_1: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 3: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_quadsect_2: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = corner_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(0);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 3: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    case (int)RefinementPattern::rp_regular: {
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = midpoint_coords.col(1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 3: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = corner_coords.col(2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      }
      break;
    }
    default: {
      LF_VERIFY_MSG(false,"No valid refinement pattern for a triangle");
      break;
    }
    } //end switch ref_pattern
    
  }

}  // namespace lf::geometry
