// This file contains Doxygen snippets for the quick reference fe_space document
// It defines a function to hold all code snippets and includes necessary
// imports.

#include <iostream>
#include <memory>

#include "lf/fe/fe.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/uscalfe/uscalfe.h"

void qr_fe_space_snippets() {
  {
    //! [fe_space Hierarchic Scalar Finite Element Space]
    // generate a simple test mesh
    std::shared_ptr<lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

    // define a polynomial degree function
    auto degree = [](const lf::mesh::Entity &entity) -> unsigned {
      return 2;  // all entities have degree 2 -> equivalent to O2 Lagrangian
                 // FESpace
    };

    // init HierarchicScalarFESpace on mesh_p
    lf::fe::HierarchicScalarFESpace<double> fe_space(mesh_p, degree);
    //! [fe_space Hierarchic Scalar Finite Element Space]
  }
  {
    //! [fe_space Lagrangian Finite Element Spaces]
    // generate a simple test mesh
    std::shared_ptr<lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(1);

    // init O1 Lagrangian finite element space on mesh_p
    lf::uscalfe::FeSpaceLagrangeO1<double> fe_space(mesh_p);
    //! [fe_space Lagrangian Finite Element Spaces]
  }
}
