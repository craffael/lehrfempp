/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 29.04.20
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
  auto mesh_file = (here.parent_path() /"/meshes/earth.msh").string();
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  /* Finite Element Space */
  auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  /* Dofhandler */
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  /* Initial Population density */
  Eigen::VectorXd u0(N_dofs); u0.setZero(); 
  u0(277) = 80;
  
  std::cout << "N_dofs :" << N_dofs << std::endl; 
  
  /* Diffusion Coefficient 
   * APPROACH 1: constant diffusion coefficient.
   *
   * auto c = [] (Eigen::Vector2d x) -> double { return 86.0;};
   *
   * APPROACH 2: Take topography into account: Mountain chains impede the dispersal of the population.
   */
  Eigen::Vector2d Himalaya(-280, 29);
  Eigen::Vector2d Alps(-350, 46);
  Eigen::Vector2d Karakoram(-285, 36);
  Eigen::Vector2d Hindukush(-288, 36);
  Eigen::Vector2d RockyMountains(-109, 44);
  Eigen::Vector2d Ural(-300, 60);
  Eigen::Vector2d Andes(-66, -21); 
  
  auto c = [&Himalaya, &Alps, &Karakoram, &Hindukush, &RockyMountains, &Ural, &Andes] (Eigen::Vector2d x) -> double {
	double diffCoeff = 86.0;
    
	if((x - Himalaya).norm() <= 3) {
	  diffCoeff = 5.0;
	  std::cout << "Take Himalaya into account." << std::endl;
	}

	if((x - Karakoram).norm() <= 3) {
      diffCoeff = 6.0;
      std::cout << "Take Karakoram into account." << std::endl;
    }
	
	if((x - Hindukush).norm() <= 2) {
      diffCoeff = 6.0;
      std::cout << "Take Hindukush into account." << std::endl;
    }

	if((x - Alps).norm() <=  3) {
      diffCoeff = 5.0;
      std::cout << "Take Alps into account." << std::endl;
    }

	if((x - Ural).norm() <= 2) {
      diffCoeff = 8.0;
      std::cout << "Take Ural into account." << std::endl;
    }

	if((x - RockyMountains).norm() <= 2) {
      diffCoeff = 10.0;
      std::cout << "Take Himalaya into account." << std::endl;
    }

    if((x - Andes).norm() <= 3) {
      diffCoeff = 9.0;
      std::cout << "Take Himalaya into account." << std::endl;
    }    

    return diffCoeff;
  };

  /* Growth Factor */
  double lambda = 1.94;

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

  /* Discern local boundary nodes */
  std::vector<Eigen::MatrixXd> locBoundary;
  for(int i = 0; i < N_dofs; ++i) {
    auto coords = lf::geometry::Corners(*(dofh.Entity(i).Geometry()));
	if(boundary_nodes(i)) {
	  locBoundary.push_back(coords);
	}
  }
  
  /* This h is to be used as a function handle for the gain term.
   * Use it in the MassEdgeMatrix Provider.
   */
  auto h = [fe_space, mesh_p, edge_pred] (Eigen::Vector2d x) -> double {
    double res = 0.0;
    /* Decaying function handle depending on x and y. */
	auto g = [x] (Eigen::Vector2d y) -> double {
	  double tmp_res = 0.0;
	  if( 5 <= (x-y).norm() &&  (x-y).norm() <= 30 ) {
		tmp_res = (1.0 / (1.0 + (x-y).squaredNorm()));
	  }
	  return tmp_res;
	};
    
	auto fe = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
    res = localQuadFunction(*mesh_p,
				{{lf::base::RefEl::kSegment(), lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2*fe->Degree())},
				{lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()}}, g, 1, edge_pred);
    
	std::cout << "res" << res << std::endl;
	return res;

  };
  
  /* In what follows, the loss term is assembled. */
  Eigen::MatrixXd L(N_dofs, N_dofs); 
  for(int j = 0; j < N_dofs; j++) {
	if(boundary_nodes(j)) {
	  auto L_j = [fe_space, &dofh, N_dofs, edge_pred, j] (Eigen::Vector2d x) -> double {
	  
	    auto g_j = [x] (Eigen::Vector2d y) -> double {
	      double tmp_res = 0.0;
		  if( 5 <= (x-y).norm() &&  (x-y).norm() <= 30 ) {
		    tmp_res = (1.0 / (1.0 + (x-y).squaredNorm()));
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

  /* Strang Splitting Method
   * SDIRK-2 evolution of linear parabolic term 
   * Exact Evolution for nonlinear reaction term
   */
  
  /* Total number of timesteps */ 
  unsigned int m = 121;
  double T = 1.; 			// the timestepsize tau will equal T/m = 0.0083
  
  /* First we assemble the carrying capacity maps */
  Eigen::MatrixXd car_cap(N_dofs, 22);

  Eigen::VectorXd K(N_dofs); K.setZero();
  K = 0.2 * Eigen::VectorXd::Ones(N_dofs);
  auto c_cap = [] (Eigen::Vector2d x) -> double { return 80.0;};
  auto h_cap = [] (Eigen::Vector2d x) -> double { return 1.0;};
  Eigen::MatrixXd L_cap(N_dofs, N_dofs); L_cap = Eigen::MatrixXd::Zero(N_dofs, N_dofs);
  /* NOTE: c = 80, lambda = 0, L = 0, h = 1, K does not matter*/
  StrangSplit DiffusionCapacity(fe_space, T, m, 0.0, c_cap, h_cap, L_cap);
  
  Eigen::VectorXd k0(N_dofs);
 
  /* t = 200 kya - 150 kya */
  k0.setZero();
  k0(277) = 80;
  car_cap.col(0) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(1) = DiffusionCapacity.Evolution(K, car_cap.col(0));
  car_cap.col(2) = DiffusionCapacity.Evolution(K, car_cap.col(1));
  car_cap.col(2) *= 10;

  std::cout << "Carrying Capacity 200kya - 150kya!" << std::endl;

  /* t = 150 kya - 130 kya */
  std::cout << "Carrying Capacity 150kya - 130kya!" << std::endl;

  /* t = 130 kya - 100 kya */
  k0.setZero();
  k0(334) = 80;
  car_cap.col(3) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(4) = DiffusionCapacity.Evolution(K, car_cap.col(3));
  car_cap.col(5) = DiffusionCapacity.Evolution(K, car_cap.col(4));
  car_cap.col(5) *= 10;

  std::cout << "Carrying Capacity 130kya - 100kya!" << std::endl;

  /* t = 100 kya - 70 kya */
  std::cout << "Carrying Capacity 100kya - 70kya!" << std::endl;

  /* t = 70 kya - 65 kya */
  k0.setZero();
  k0(253) = 80;
  car_cap.col(6) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(7) = DiffusionCapacity.Evolution(K, car_cap.col(6));
  car_cap.col(7) *= 15;

  std::cout << "Carrying Capacity 70kya - 65kya!" << std::endl;

  /* t = 65 kya - 50 kya */
  k0.setZero();
  k0(1775) = 80;
  k0(222) = 80;
  k0(181) = 80;
  car_cap.col(8) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(9) = DiffusionCapacity.Evolution(K, car_cap.col(8));
  car_cap.col(9) *= 8;
  
  std::cout << "Carrying Capacity 65kya - 50kya!" << std::endl;
 
  /* t = 50 kya - 45 kya */
  k0.setZero();
  k0(400) = 200;
  car_cap.col(10) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(10) *= 2;

  std::cout << "Carrying Capacity 50kya - 45kya!" << std::endl;

  /* t = 45 kya - 25 kya */
  k0.setZero();
  k0(492) = 80;
  k0(114) = 80;
  car_cap.col(11) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(12) = DiffusionCapacity.Evolution(K, car_cap.col(11));
  car_cap.col(12) *= 8;

  std::cout << "Carrying Capacity 45kya - 25kya!" << std::endl;

  /* t = 25 kya - 15 kya */
  k0.setZero();
  k0(1342) = 80;
  car_cap.col(13) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(14) = DiffusionCapacity.Evolution(K, car_cap.col(13));
  car_cap.col(14) *= 25;
 
  std::cout << "Carrying Capacity 25kya - 15kya!" << std::endl;

  /* t = 15 kya - 0 kya */
  k0.setZero();
  k0(1257) = 80;
  k0(1100) = 80;
  car_cap.col(15) = DiffusionCapacity.Evolution(K, k0);
  car_cap.col(15) *= 4;

  std::cout << "Carrying Capacity 15kya - 0kya!" << std::endl;
  std::cout << "Capacity maps are assembled!" << std::endl;

  /* Now we may compute the solution */
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c, h, L);
  
  Eigen::MatrixXd sol(N_dofs, 40);
  Eigen::VectorXd cap(N_dofs); cap.setZero();
  /* t = 200 kya - 150 kya */
  cap = car_cap.col(2);
  sol.col(0) = StrangSplitter.Evolution(cap, u0);
  
  std::cout << "Solution 200kya - 150kya!" << std::endl;
  
  /* t = 150 kya - 130 kya */
  sol.col(1) = StrangSplitter.Evolution(cap, sol.col(0));

  std::cout << "Solution 150kya - 130kya!" << std::endl;

  /* t = 130 kya - 100 kya */
  cap = cap + car_cap.col(5);
  sol.col(2) = StrangSplitter.Evolution(cap, sol.col(1));
  
  std::cout << "Solution 130kya - 100kya!" << std::endl;

  /* t = 100 kya - 70 kya */
  sol.col(3) = StrangSplitter.Evolution(cap, sol.col(2));
  
  std::cout << "Solution 100kya - 70kya!" << std::endl;

  /* t = 70 kya - 65 kya */
  cap = cap + car_cap.col(7); 
  sol.col(4) = StrangSplitter.Evolution(cap, sol.col(3));
  sol.col(5) = StrangSplitter.Evolution(cap, sol.col(4));

  std::cout << "Solution 70kya - 65kya!" << std::endl;

  /* t = 65 kya - 50 kya */
  cap = cap + car_cap.col(9);
  sol.col(6) = StrangSplitter.Evolution(cap, sol.col(5));
  sol.col(7) = StrangSplitter.Evolution(cap, sol.col(6));
  
  std::cout << "Solution 65kya - 50kya!" << std::endl;

  /* t = 50 kya - 45 kya */
  cap = cap + car_cap.col(10);
  sol.col(8) = StrangSplitter.Evolution(cap, sol.col(7));
  
  std::cout << "Solution 50kya - 45kya!" << std::endl;

  /* t = 45 kya - 25 kya */
  cap = cap + car_cap.col(12);
  sol.col(9) = StrangSplitter.Evolution(cap, sol.col(8));
  sol.col(10) = StrangSplitter.Evolution(cap, sol.col(9));
  sol.col(11) = StrangSplitter.Evolution(cap, sol.col(10));
  sol.col(12) = StrangSplitter.Evolution(cap, sol.col(11));

  std::cout << "Solution 45kya - 25kya!" << std::endl;

  /* t = 25 kya - 15 kya */
  cap = cap + car_cap.col(14);
  sol.col(13) = StrangSplitter.Evolution(cap, sol.col(12));
  sol.col(14) = StrangSplitter.Evolution(cap, sol.col(13));
  
  std::cout << "Solution 25kya - 15kya!" << std::endl;

  /* t = 15 kya - 0 kya */
  cap = cap + car_cap.col(15);
  sol.col(15) = StrangSplitter.Evolution(cap, sol.col(14));
  sol.col(16) = StrangSplitter.Evolution(cap, sol.col(15));
  sol.col(17) = StrangSplitter.Evolution(cap, sol.col(16));
  sol.col(18) = StrangSplitter.Evolution(cap, sol.col(17));
  sol.col(19) = StrangSplitter.Evolution(cap, sol.col(18));
  
  std::cout << "Solution 15kya - 0kya!" << std::endl;

  /* Use VTK-Writer for Visualization of solution */
  lf::io::VtkWriter vtk_writer_cap(mesh_p, "cap.vtk");
  auto nodal_data_cap = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data_cap->operator()(dofh.Entity(global_idx)) = cap[global_idx];
  }
  vtk_writer_cap.WritePointData("cap", *nodal_data_cap);

  for(int k = 1; k < 21; k++) {

    std::stringstream filename;
    filename << "sol" << k << ".vtk";

    lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data->operator()(dofh.Entity(global_idx)) = sol.col(k-1)[global_idx];
    }
    
	vtk_writer.WritePointData("sol", *nodal_data);
  }

  for(int k = 1; k < 17; k++) {

    std::stringstream filename_cap;
    filename_cap << "cap" << k << ".vtk";

    lf::io::VtkWriter vtk_writer_caps(mesh_p, filename_cap.str());
    auto nodal_data_caps = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data_caps->operator()(dofh.Entity(global_idx)) = car_cap.col(k-1)[global_idx];
    }
    
	vtk_writer_caps.WritePointData("caps", *nodal_data_caps);
  }


  return 0;
}
