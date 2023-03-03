/**
 * @file
 * @headerfile projects/dpg/dpg_tools.h
 * @brief implementation of dpg tools
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */
#include "dpg_tools.h"

namespace projects::dpg {

std::vector<lf::quad::QuadRule> BoundaryQuadRule(lf::base::RefEl ref_el,
                                                 const lf::quad::QuadRule& qr) {
  // construct the geometry of the reference element
  std::shared_ptr<lf::geometry::Geometry> geo_ptr;
  switch (ref_el) {
    case (lf::base::RefEl::kTria()):
      geo_ptr = std::make_shared<lf::geometry::TriaO1>(
          lf::base::RefEl::kTria().NodeCoords());
      break;
    case (lf::base::RefEl::kQuad()):
      geo_ptr = std::make_shared<lf::geometry::QuadO1>(
          lf::base::RefEl::kQuad().NodeCoords());
      break;
    default:
      LF_ASSERT_MSG(false, ref_el << "unsupported reference element.");
  }

  // query the number of edges
  lf::base::size_type numSegments = ref_el.NumSubEntities(1);
  std::vector<lf::quad::QuadRule> BoundaryQuadRules(numSegments);

  // iterate over edges
  for (int segment = 0; segment < numSegments; ++segment) {
    // query edge geometry and transform weights
    auto edge_ptr = geo_ptr->SubGeometry(1, segment);
    Eigen::MatrixXd Points = edge_ptr->Global(qr.Points());
    Eigen::VectorXd Weights =
        qr.Weights() *
        edge_ptr->IntegrationElement((Eigen::VectorXd(1) << 0.5).finished());
    lf::quad::QuadRule bqr(ref_el, Points, Weights, 0);
    BoundaryQuadRules[segment] = bqr;
  }

  return BoundaryQuadRules;
}

Eigen::MatrixXd OuterNormals(const lf::geometry::Geometry& geometry) {
  // check, that the reference element is supported:
  lf::base::RefEl refEl = geometry.RefEl();
  LF_ASSERT_MSG(
      refEl == lf::base::RefEl::kTria() || refEl == lf::base::RefEl::kQuad(),
      "invalid reference element: " << refEl);

  // retrive the corners of the geometry object.
  Eigen::MatrixXd corners = lf::geometry::Corners(geometry);
  lf::base::size_type n_corners = corners.cols();
  lf::base::size_type n_edges = n_corners;

  // determine, if the numbering of corners is
  // counterclockwise or clockwise. This is needed to conclude
  // in which direction the edge vectors have to be turned to
  // retrive the OUTER normal.
  int orientation;
  if (refEl == lf::base::RefEl::kTria()) {
    // check sign of the JacobianDeterminant to deduce
    // if the orientation of the numbering has changed.
    // since the Jacobian is constant, it can be evaluated in any
    // point of the reference triangle.
    Eigen::MatrixXd point = Eigen::MatrixXd::Zero(2, 1);
    Eigen::MatrixXd Jacobian = geometry.Jacobian(point);
    double JacobianDeterminant = Jacobian.determinant();
    orientation = JacobianDeterminant > 0 ? 1 : -1;

  } else {
    // in the case of a quadrilateral element it can
    // happen, that the resulting global quadrialteral
    // is non convex, then evaluation of the (non constant)
    // Jacobian determinant in  one point is no longer enough to determine the
    // orientation. instead the Jacobian is evaluated in all corners of the
    // reference element.
    Eigen::MatrixXd Jacobians = geometry.Jacobian(refEl.NodeCoords());

    // Count the number of Jacobian evaluations with a positive determinant.
    int positiveJacobianDeterminantCount = 0;
    for (int i = 0; i < 4; ++i) {
      Eigen::MatrixXd Jacobian = Jacobians.block(0, 2 * static_cast<Eigen::Index>(i), 2, 2);
      double JacobianDeterminant = Jacobian.determinant();
      if (JacobianDeterminant > 0) {
        positiveJacobianDeterminantCount++;
      }
    }
    // Jacobian determinant is positive in all points:
    //--> convex quadrilateral, counterclockwise orientation.
    // Jacobian determinant is positive in all but one point:
    //--> concave quadrilateral, counterclickwise orientation.
    // Jacobian determinant negative in all points:
    //--> convex quadrilateral, clockwise orientation.
    // Jacobian determinant positive in only one point:
    //--> concave quadrilateral, clockwise orientationi.
    orientation = positiveJacobianDeterminantCount >= 3 ? 1 : -1;
  }

  // calculate  "edge" direction vectors:
  Eigen::MatrixXd edges(2, n_edges);
  for (int i = 0; i < n_edges; ++i) {
    edges.col(i) = corners.col(refEl.SubSubEntity2SubEntity(1, i, 1, 1)) -
                   corners.col(refEl.SubSubEntity2SubEntity(1, i, 1, 0));
  }

  // calculate normal vectors and normalize them.
  // take into account the possible clockwise orientation of numbering on
  // the triangle, in which case the normal vector has to be flipped.
  Eigen::MatrixXd normals(2, n_edges);
  for (int i = 0; i < n_edges; ++i) {
    Eigen::Vector2d normal;
    normal(0) = edges(1, i) * orientation;
    normal(1) = -edges(0, i) * orientation;
    normal.normalize();
    normals.col(i) = normal;
  }
  return normals;
}

// initializes the maximum element data set:*/
void PrescribedSignProvider::init() {
  LF_ASSERT_MSG((mesh_ptr_ != nullptr), "invalid mesh (nullptr)");

  lf::mesh::utils::CodimMeshDataSet<int>& maxElement = *maxElement_ptr_;
  // find the maximum index of cells to which any given edge belongs.
  for (const lf::mesh::Entity* const e : mesh_ptr_->Entities(0)) {
    auto entity_idx = mesh_ptr_->Index(*e);
    for (const lf::mesh::Entity* const subent : e->SubEntities(1)) {
      LF_ASSERT_MSG(maxElement.DefinedOn(*subent),
                    "maxElement_ not defined on subentity" << *subent);
      maxElement(*subent) = std::max<int>(maxElement(*subent), lf::base::narrow<int>(entity_idx));
    }
  }
}

// evaluate the prescribed sign using the described maximum rule.
int PrescribedSignProvider::PrescribedSign(const lf::mesh::Entity& element,
                                           const lf::mesh::Entity& edge) const {
  lf::mesh::utils::CodimMeshDataSet<int>& maxElement = *maxElement_ptr_;
  LF_ASSERT_MSG(maxElement.DefinedOn(edge),
                "maxElement_ not defined on edge " << edge);
  return mesh_ptr_->Index(element) == maxElement(edge) ? 1 : -1;
}

}  // namespace projects::dpg
