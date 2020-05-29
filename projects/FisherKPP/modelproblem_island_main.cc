/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 22.04.20
 *  @copyright ETH Zurich
 */

#include "strangsplitting.cc"

#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

#include <lf/io/io.h>

using namespace FisherKPP;

int main(int /*argc*/, char ** /*argv*/){

  /* Obtain mesh */
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  boost::filesystem::path here = __FILE__;
  auto mesh_file = (here.parent_path() /"/meshes/Island.msh").string();
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  /* Finite Element Space */
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  /* Dofhandler */
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  /* Initial Population density */
  Eigen::VectorXd u0(N_dofs); u0.setZero(); 
  u0(51) = 0.3;
  std::cout << "N_dofs :" << N_dofs << std::endl; 
  
  /* Diffusion Coefficient */
  Eigen::Vector2d init(-330, -14);
  // TEST: non constant diffusion coefficient: take into account topology.
  auto c = [&init] (Eigen::Vector2d x) -> double { 
	double diffcoeff = 1.2;

    if((x - init).norm() <= 1) {
	  diffcoeff = 0.2;
	  std::cout << "Take topography into account." << std::endl;
	}

	return diffcoeff;
  };
  
  /* Growth Factor */
  double lambda = 2.1;

  /* Non Local Boundary Conditions */

  /* This predicate returns true for nodes on the boundary */
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  auto boundary_nodes = [&bd_flags, &dofh] (unsigned int idx) -> bool {
      return bd_flags(dofh.Entity(idx));
  };
  
  /* This predicate returns true for edges on the boundary */
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  auto edge_pred = [&bd_flags_edge] (const lf::mesh::Entity &edge) -> bool {
	  return bd_flags_edge(edge);
  };
  
  
  /* This h is to be used as a function handle for the gain term.
   * Use it in the MassEdgeMatrix Provider.
   */
  auto h = [fe_space, mesh_p, edge_pred] (Eigen::Vector2d x) -> double {
    double res = 0.0;
    /* Decaying function handle depending on x and y.
     * g has a non negative contribution to the integral value
     * for all nodes which are not too far away from x.
     */
	auto g = [x] (Eigen::Vector2d y) -> double {
	  double tmp_res = 0.0;
	  // NORM 
	  if((x-y).norm() >= 15 &&(x-y).norm() <= 35 ) {
		tmp_res = ( 1.0 / (1.0 + (x-y).squaredNorm()) );
	 }
	  
      return tmp_res;
	
	};
    
	auto fe = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
    res = localQuadFunction(*mesh_p,
				{{lf::base::RefEl::kSegment(), lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2*fe->Degree())},
				{lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()}}, g, 1, edge_pred);
    std::cout << "rse " << res << std::endl;
	return res;

  };
  /* In what follows, the loss term is assembled. */
  Eigen::MatrixXd L(N_dofs, N_dofs); 
  for(int j = 0; j < N_dofs; j++) {
    
	if(boundary_nodes(j)) {

	auto L_j = [fe_space, &dofh, N_dofs, edge_pred, j] (Eigen::Vector2d x) -> double {
	  
 	  auto g_j = [x] (Eigen::Vector2d y) -> double {
        double tmp_res = 0.0;
		if((x-y).norm() >= 15 &&(x-y).norm() <= 35) {
		  tmp_res = (1.0 / (1.0 + (x-y).squaredNorm()) );
		}
		
		return tmp_res;
	   
	 };
	 
	  Eigen::VectorXd L1(N_dofs); L1.setZero();

	  lf::mesh::utils::MeshFunctionGlobal mf_g{g_j};
	  lf::uscalfe::ScalarLoadEdgeVectorProvider<double, decltype(mf_g), decltype(edge_pred)>  edgeVec_y(fe_space, mf_g, edge_pred);
	  lf::assemble::AssembleVectorLocally(1, dofh, edgeVec_y, L1);
	 
	  return L1(j);
  
    };
    
	Eigen::VectorXd L2(N_dofs); L2.setZero();

    lf::mesh::utils::MeshFunctionGlobal mf_L{L_j};
  	lf::uscalfe::ScalarLoadEdgeVectorProvider<double, decltype(mf_L), decltype(edge_pred)>  edgeVec_x(fe_space, mf_L, edge_pred);
	lf::assemble::AssembleVectorLocally(1, dofh, edgeVec_x, L2);
    
	L.col(j) = L2;
	} 

	else {
	  L.col(j) = Eigen::VectorXd::Zero(N_dofs);
	}

  }
  std::cout << "L " << L.norm() << std::endl;
  
  /* Strang Splitting Method
   * SDIRK-2 evolution of linear parabolic term 
   * Exact Evolution for nonlinear reaction term
   */
  
  /* Total number of timesteps */ 
  unsigned int m = 121;
  double T = 1.; 			// the timestepsize tau will equal T/m = 0.0083
  

  /* Now we may compute the solution */
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c, h, L);
  
  Eigen::MatrixXd sol(N_dofs, 12);
  Eigen::VectorXd cap(N_dofs); cap.setZero();
  cap = 0.9 * Eigen::VectorXd::Ones(N_dofs);
  
  sol.col(0) = StrangSplitter.Evolution(cap, u0);
  for(int i = 1; i < 10; i++) {
    sol.col(i) = StrangSplitter.Evolution(cap, sol.col(i-1));
  }

  /* Use VTK-Writer for Visualization of solution */
  for(int k = 1; k < 11; k++) {

    std::stringstream filename;
    filename << "sol" << k << ".vtk";

    lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data->operator()(dofh.Entity(global_idx)) = sol.col(k-1)[global_idx];
    }
    
	vtk_writer.WritePointData("sol", *nodal_data);
  }

  return 0;

}
