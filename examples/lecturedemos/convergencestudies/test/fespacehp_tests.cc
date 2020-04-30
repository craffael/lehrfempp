#define _USE_MATH_DEFINES

#include <gtest/gtest.h>

#include "../fespacehp.h"

#include <cmath>
#include <array>
#include <tuple>
#include <string>
#include <memory>
#include <lf/base/base.h>
#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/refinement/refinement.h>
#include <lf/refinement/mesh_function_transfer.h>
#include <lf/io/io.h>


/**
 * @brief Builds a mesh on [0, 1]^2 from two triangles
 * @returns A shared pointer to a mesh
 */
std::shared_ptr<lf::mesh::Mesh> getSquareDomainTriangles() {
  lf::mesh::hybrid2d::MeshFactory factory(2);
  // Add the vertices
  std::vector<lf::mesh::MeshFactory::size_type> vertices;
  Eigen::Vector2d vertex_coord;
  vertex_coord << 0, 0;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 1, 0;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 1, 1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 0, 1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  // Add the triangles
  Eigen::Matrix<double, Eigen::Dynamic, 3> coords(2, 3);
  lf::mesh::MeshFactory::size_type nodes[3];
  coords << 0, 1, 0, 0, 0, 1;
  nodes[0] = vertices[0];
  nodes[1] = vertices[1];
  nodes[2] = vertices[3];
  auto geom_tria1 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria1));
  coords << 1, 1, 0,
	    0, 1, 1;
  nodes[0] = vertices[1];
  nodes[1] = vertices[2];
  nodes[2] = vertices[3];
  auto geom_tria2 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria2));
  // Build the mesh
  return factory.Build();
}


/**
 * @brief Builds a mesh on [0, 1]^2 from one quad
 * @returns A shared pointer to a mesh
 */
std::shared_ptr<lf::mesh::Mesh> getSquareDomainQuad() {
  lf::mesh::hybrid2d::MeshFactory factory(2);
  // Add the vertices
  std::vector<lf::mesh::MeshFactory::size_type> vertices;
  Eigen::Vector2d vertex_coord;
  vertex_coord << 0, 0;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 1, 0;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 1, 1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 0, 1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  // Add the quad
  Eigen::Matrix<double, Eigen::Dynamic, 4> coords(2, 4);
  lf::mesh::MeshFactory::size_type nodes[4];
  coords << 0, 1, 1, 0,
	    0, 0, 1, 1;
  nodes[0] = vertices[0];
  nodes[1] = vertices[1];
  nodes[2] = vertices[2];
  nodes[3] = vertices[3];
  auto geom_quad = std::make_unique<lf::geometry::QuadO1>(coords);
  factory.AddEntity(lf::base::RefEl::kQuad(), nodes, std::move(geom_quad));
  // Build the mesh
  return factory.Build();
}


TEST(fespacelagrangeon, bilinear_laplacian) {
    // Evaluate the bilinear form at this function to test the convergence
    const auto w = [](const Eigen::Vector2d& x) {
	return std::sin(M_PI*x[0]) * std::sin(M_PI*x[1]);
    };
    lf::mesh::utils::MeshFunctionGlobal mf_w(w);
    // The exact result of b(w, w)
    const double bilinear_exact = M_PI*M_PI/2;
    // Generate the meshes on which to test the bilinear form
    std::array<std::tuple<std::string, std::shared_ptr<lf::mesh::Mesh>>, 2> meshes{
	std::tuple<std::string, std::shared_ptr<lf::mesh::Mesh>>{"Triangles", getSquareDomainTriangles()},
	std::tuple<std::string, std::shared_ptr<lf::mesh::Mesh>>{"Quad", getSquareDomainQuad()}
    };

    // Test the convergence for each mesh
    const unsigned max_p = 20;
    for (const auto& [name, mesh] : meshes) {
	Eigen::VectorXd error(max_p);
	for (unsigned p = 1 ; p <= max_p ; ++p) {
	    // Construct the FE space
	    const auto fe_space = std::make_shared<FeSpaceHP<double>>(mesh, p);
	    const auto& dofh = fe_space->LocGlobMap();
	    // Assemble the bilinear form
	    const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
	    const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(0);
	    lf::uscalfe::ReactionDiffusionElementMatrixProvider elem_mat_prov(fe_space, mf_alpha, mf_gamma);
	    lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
	    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_prov, A_COO);
	    Eigen::SparseMatrix<double> A = A_COO.makeSparse();
	    // Interpolate w on the FE space
	    const auto w_dofs = lf::uscalfe::NodalProjection(*fe_space, mf_w);
	    // Compute the discretized bilinear form
	    const double bilinear_approx = w_dofs.transpose() * A * w_dofs;
	    // Store the error
	    error[p-1] = bilinear_exact - bilinear_approx;
	    // Store the interpolation to vtk for debugging purposes
	    /*
	    const lf::uscalfe::MeshFunctionFE<double, double> mf_w_interp(fe_space, w_dofs);
	    const auto mh = lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, 6);
	    const lf::refinement::MeshFunctionTransfer mf_w_fine(*mh, mf_w_interp, 0, 6);
	    lf::io::VtkWriter writer(mh->getMesh(6), name + "_" + std::to_string(p) + "_w_interp.vtk");
	    writer.WritePointData("exact", mf_w);
	    writer.WritePointData("interp", mf_w_fine);
	    */
	}
	EXPECT_TRUE(false) << "Mesh " << name << "\n" << error.transpose() << std::endl;
    }
}
