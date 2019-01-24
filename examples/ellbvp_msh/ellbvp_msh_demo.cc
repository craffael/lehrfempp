/**
 * @file
 * @brief Solution of general second-order elliptic boundary value problem with
 * linear finite elements from a Gmsh generated mesh
 * @author Simon Meierhans
 * @date   January 2019
 * @copyright MIT License
 */

#include <fstream>
#include <iomanip>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/filesystem.hpp>


int main(){
  // abbreviations for types
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;

  //find path to mesh
  boost::filesystem::path here = __FILE__;
  auto mesh_path = here.parent_path() / "meshes/square.msh";

  //load mesh of square computational domain
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());

  //get pointer to mesh
  auto mesh = reader.mesh();

  // Count the number of edges with von Neumann boundary condition
  auto physical_entity_nr_neu = reader.PhysicalEntityName2Nr("neu");
  int num_neumann = 0;
  for (auto& e : mesh->Entities(1)) {
    if (reader.IsPhysicalEntity(e, physical_entity_nr_neu)) {
      ++num_neumann;
    }
  }
  std::cout << "boundary edges with von Neumann boundary condition: " << num_neumann << "\n";

  // Count the number of edges with Dirichlet boundary condition
  auto physical_entity_nr_dir = reader.PhysicalEntityName2Nr("dir");
  int num_dirichlet = 0;
  for (auto& e : mesh->Entities(1)) {
    if (reader.IsPhysicalEntity(e, physical_entity_nr_dir)) {
      ++num_dirichlet;
    }
  }
  std::cout << "boundary edges with Dirichlet boundary condition: " << num_dirichlet << "\n";

  //set up finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  //set up dof handler
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // Dimension of finite element space`
  const size_type N_dofs(dofh.NoDofs());

  //identity mesh function for very simple problem
  auto identity = [](Eigen::Vector2d x) -> double {
    return 1.;
  };
  lf::uscalfe::MeshFunctionGlobal mf_identity{identity};

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);

  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(mf_identity), decltype(mf_identity)>
      elmat_builder(fe_space, mf_identity, mf_identity);

  // Invoke assembly on cells (co-dimension = 0)
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);

  // Right-hand side vector; has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // Initialize object taking care of local computations on all cells for the source f. The source is the identity function
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_identity)>
      elvec_builder(fe_space, mf_identity);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  //Select von Neumann edges
  neu_edge_sel

}
