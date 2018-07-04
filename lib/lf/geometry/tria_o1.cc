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
  TriaO1::ChildGeometry(int ref_pattern,int anchor,int selector) const {
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

    // Remap local indices according to anchor values
    const int mod_0 = (0+anchor)%3;
    const int mod_1 = (1+anchor)%3;
    const int mod_2 = (2+anchor)%3;
    
    // Create child geometries according to refinement patterns and selection
    switch (ref_pattern) {
    case (int)RefinementPattern::rp_nil: {
      return(nullptr);
      break;
    }
    case (int)RefinementPattern::rp_copy: {
      return std::make_unique<TriaO1>(coords_);
      break;
    }
    case (int)RefinementPattern::rp_bisect: {
      // Splitting a triangle in two by bisecting anchor edge
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(mod_0);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = corner_coords.col(mod_2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(mod_1);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = corner_coords.col(mod_2);
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
    case (int)RefinementPattern::rp_trisect: {
      // Bisect through anchor edge first and then bisect through
      // edge with the next larger index (mod 3); creates three
      // child triangles.
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(mod_0);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = corner_coords.col(mod_2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(mod_1);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(mod_2);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_2);
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
    case (int)RefinementPattern::rp_trisect_left: {
      // Bisect through anchor edge first and then bisect through
      // edge with the next smaller index (mod 3); creates three
      // child triangles.
       switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(mod_0);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(mod_1);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = corner_coords.col(mod_2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(mod_2);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_2);
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
    case (int)RefinementPattern::rp_quadsect: {
      // Bisect through the anchor edge first and then
      // through the two remaining edges; creates four child
      // triangles; every edge is split.
      switch(selector) {
      case 0: { 
	child_coords.col(0) = corner_coords.col(mod_0);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_2);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(mod_1);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(mod_2);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_1);
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 3: {
	child_coords.col(0) = corner_coords.col(mod_2);
	child_coords.col(1) = midpoint_coords.col(mod_0);
	child_coords.col(2) = midpoint_coords.col(mod_2);
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
      // Split triangle into four small  congruent triangles
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
    case (int)RefinementPattern::rp_barycentric: {
      //  Split triangle into 5 smaller triangles by connecting
      // the center of gravity with the vertices and the midpoints
      // of the edges.
      
      // Obtain coordinates of center of gravity
      Eigen::Matrix<double,2,1> ref_baryc_coords = Eigen::Vector2d({1.0/3,1.0/3});
      Eigen::MatrixXd baryc_coords = Global(ref_baryc_coords);
      switch(selector) {
      case 0: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = baryc_coords;
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 1: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(0);
	child_coords.col(2) = baryc_coords;
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 2: {
	child_coords.col(0) = corner_coords.col(1);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = baryc_coords;
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 3: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(1);
	child_coords.col(2) = baryc_coords;
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 4: {
	child_coords.col(0) = corner_coords.col(2);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = baryc_coords;
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      case 5: {
	child_coords.col(0) = corner_coords.col(0);
	child_coords.col(1) = midpoint_coords.col(2);
	child_coords.col(2) = baryc_coords;
	return std::make_unique<TriaO1>(child_coords);
	break;
      }
      default: {
	LF_VERIFY_MSG(false,"Selector " << selector
		      << " invalid for pattern " << ref_pattern);
	break;
      }
      } // end switch selector
      break;
    }
    default: {
      LF_VERIFY_MSG(false,"No valid refinement pattern for a triangle");
      break;
    }
    } //end switch ref_pattern
    
  }

}  // namespace lf::geometry
