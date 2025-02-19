// This file contains Doxygen snippets for the quick reference geometry document
// It defines a function to hold all code snippets and includes necessary
// imports.

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>

#include <iostream>

void qr_geometry_snippets() {
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);
  unsigned codim = 0;

  {
    //! [geometry Geometry Interface]
    for (const lf::mesh::Entity* entity : mesh_p->Entities(codim)) {
      // Get the geometry of an entity
      const lf::geometry::Geometry* geometry = entity->Geometry();
    }
    //! [geometry Geometry Interface]
  }

  auto entity = mesh_p->EntityByIndex(0, 0);

  {
    //! [geometry Utility Functions]
    const lf::geometry::Geometry* geometry = entity->Geometry();

    // Get the volume of an entity (for edges length, for triangles area, etc.)
    double v = lf::geometry::Volume(*geometry);

    // Get the corners of an entity as a dxn matrix
    // where d is the dimension of physical space, and n the number of corners.
    // (for edges 2 corners, for triangles 3 corners, etc.)
    Eigen::MatrixXd corners = lf::geometry::Corners(*geometry);
    //! [geometry Utility Functions]

    //! [geometry Utility Functions 0]
    // Test if geometry is a non-degenerate bilinear quadrilateral (corners
    // matrix with 4 cols)
    bool is_not_deg_q = lf::geometry::assertNonDegenerateQuad(corners);

    // Test if geometry is a non-degenerate triangle (corners matrix with 3
    // cols)
    bool is_not_deg_t = lf::geometry::assertNonDegenerateTriangle(corners);
    //! [geometry Utility Functions 0]
  }

  auto geometry = entity->Geometry();

  {
    //! [geometry Geometry Methods]
    // Dimension of the local coordinate system
    unsigned dim_local = geometry->DimLocal();

    // Dimension of the physical coordinate system
    unsigned dim_global = geometry->DimGlobal();

    // Get the reference element for the geometry
    lf::base::RefEl ref_el = geometry->RefEl();

    // Check if mapping is affine
    bool is_affine = geometry->isAffine();
    //! [geometry Geometry Methods]
  }

  {
    //! [geometry Global]
    // Define two point on the reference element (local).
    Eigen::MatrixXd local_points(2, 2);
    local_points << 0.1, 0.5, 0.6, 0.2;

    // Map the points from the local into the global space.
    Eigen::MatrixXd global = geometry->Global(local_points);
    //! [geometry Global]

    //! [geometry IntegrationElement]
    // Compute the integration element for each point in points
    Eigen::VectorXd integration_element =
        geometry->IntegrationElement(local_points);
    //! [geometry IntegrationElement]

    //! [geometry Jacobian]
    // Compute the Jacobian matrix for each point in points
    Eigen::MatrixXd jacobian = geometry->Jacobian(local_points);
    //! [geometry Jacobian]

    unsigned n = 1;
    //! [geometry Jacobian 0]
    unsigned dim_global = geometry->DimGlobal();
    unsigned dim_local = geometry->DimLocal();

    // Access the Jacobian matrix for the n-th point (starting at 0)
    Eigen::MatrixXd jacobian_n =
        jacobian.block(0, n * dim_local, dim_global, dim_local);
    //! [geometry Jacobian 0]

    //! [geometry Jacobian 1]
    // Access 2x2 Jacobian matrix starting at row 0 and column 2
    Eigen::MatrixXd jacobian_2 = jacobian.block(0, 2, 2, 2);
    //! [geometry Jacobian 1]

    //! [geometry JacobianInverseGramian]
    // Compute the inverse of the Jacobian matrix for each point in points
    Eigen::MatrixXd jacobian_inv =
        geometry->JacobianInverseGramian(local_points);

    // Access the JacobianInverseGramian matrix for the n-th point (starting
    // at 0)
    Eigen::MatrixXd jacobian_inv_n =
        jacobian_inv.block(0, n * dim_local, dim_global, dim_local);
    //! [geometry JacobianInverseGramian]
  }
}
