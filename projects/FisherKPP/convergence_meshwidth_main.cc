/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 11.05.20
 *  @copyright ETH Zurich
 */

#include "norms.cc"
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
 
  Eigen::VectorXd numDofs(5); numDofs.setZero();
  Eigen::VectorXd meshsizes(5); meshsizes.setZero();
  
  Eigen::VectorXd eL2(5); eL2.setZero();
  Eigen::VectorXd u_ref(6474); 
  
  Eigen::VectorXd mu5(6474); mu5.setZero();
  Eigen::VectorXd mu4(1670); mu4.setZero();
  Eigen::VectorXd mu3(444); mu3.setZero();
  Eigen::VectorXd mu2(125); mu2.setZero();
  Eigen::VectorXd mu1(39); mu1.setZero();
  
  Eigen::VectorXd mu_sub1(39); mu_sub1.setZero();
  Eigen::VectorXd mu_sub2(125); mu_sub2.setZero();
  Eigen::VectorXd mu_sub3(444); mu_sub3.setZero();
  Eigen::VectorXd mu_sub4(1670); mu_sub4.setZero();
  Eigen::VectorXd mu_sub5(6474); mu_sub5.setZero();

  /* Diffusion Coefficient */
  auto c = [] (Eigen::Vector2d x) -> double { return 1.2;};
  /* Growth Factor */
  double lambda = 2.1;
  
  for(int l = 1; l <= 5; l++) {
    
	/* Obtain mesh */
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    boost::filesystem::path here = __FILE__;
	auto mesh_file = (here.parent_path() / ("/meshes/test" + std::to_string(l) + ".msh")).string();
	lf::io::GmshReader reader(std::move(mesh_factory), mesh_file); 
	std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

    /* Finite Element Space */
    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    /* Dofhandler */
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
    numDofs(l-1) = N_dofs;
	meshsizes(l-1) = getMeshSize(mesh_p);

	/* Initial Population density */
    Eigen::VectorXd u0(N_dofs); u0.setZero(); 
    
	if(l == 1) {
	  u0(13) = 0.3;
	} else if(l == 2) {
	  u0(26) = 0.3;
	} else if(l == 3) {
	  u0(52) = 0.3;
	} else if(l == 4) {
      u0(104) = 0.3;
    } else if(l == 5) {
      u0(208) = 0.3;
    }

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
	  
	  auto g = [x] (Eigen::Vector2d y) -> double {
	    double tmp_res = 0.0;
	    if((x-y).norm() >= 15 &&(x-y).norm() <= 35 ) {
		  tmp_res = ( 1.0 / (1.0 + (x-y).squaredNorm()) );
	    }
        return tmp_res;
	  };
    
	  auto fe = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
      res = localQuadFunction(*mesh_p,
				{{lf::base::RefEl::kSegment(), lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2*fe->Degree())},
				{lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()}}, g, 1, edge_pred);
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
    std::cout << "L norm " << L.norm() << std::endl;
  
    /* Strang Splitting Method
     * SDIRK-2 evolution of linear parabolic term 
     * Exact Evolution for nonlinear reaction term
     */
  
    /* Total number of timesteps */ 
    unsigned int m = 100;
    double T = 1.; 			// the timestepsize tau will equal T/m = 0.01
	
    /* Now we may compute the solution */
    StrangSplit StrangSplitter(fe_space, T, m, lambda, c, h, L);

	Eigen::VectorXd cap(N_dofs); cap.setZero();
    cap = 0.8 * Eigen::VectorXd::Ones(N_dofs);
    
	Eigen::VectorXd sol(N_dofs);
    sol = StrangSplitter.Evolution(cap, u0);
    
	if(l == 5) {
	  // Reference solution: solution on finest mesh.
	  u_ref = sol;
	}
    
	if(l == 1) {
	  mu1 = sol;
    } else if(l == 2) {
	  mu2 = sol;
    } else if(l == 3) {
      mu3 = sol;
	} else if(l == 4) {
	  mu4 = sol;
	}else if(l == 5) {
	  mu5 = sol;
	}

	std::cout << mesh_file << std::endl;
	std::cout << "\t(Ndofs = " << N_dofs << "):" << std::endl;
    std::cout << std::setprecision(16);

  }
  
  /* Compare the solution to the Reference solution on the finest mesh. */
  mu_sub1 = reduce(u_ref, numDofs(0));
  mu_sub2 = reduce(u_ref, numDofs(1));
  mu_sub3 = reduce(u_ref, numDofs(2));
  mu_sub4 = reduce(u_ref, numDofs(3));
  mu_sub5 = reduce(u_ref, numDofs(4));

  eL2(0) = (mu1 - mu_sub1).lpNorm<2>();
  eL2(1) = (mu2 - mu_sub2).lpNorm<2>();
  eL2(2) = (mu3 - mu_sub3).lpNorm<2>();
  eL2(3) = (mu4 - mu_sub4).lpNorm<2>();
  eL2(4) = (mu5 - mu_sub5).lpNorm<2>();
  
  for(int l = 0; l < 5; l++) {

	std::cout << "numdofs for l = " << l << " : " << numDofs(l) << std::endl;
	std::cout << "meshsize for l = " << l << " : " << meshsizes(l) << std::endl;
    
    std::cout << "L2 error of solution with respect to reference solution on finest mesh, for l = " << l << " : " << eL2(l) << std::endl;
    
  }

  // Define output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  
  std::ofstream file;

  file.open("errorL2.csv");
  file << eL2.format(CSVFormat);
  file.close();
  
  file.open("NumDofs.csv");
  file << numDofs.format(CSVFormat);
  file.close();

  file.open("meshsizes.csv");
  file << meshsizes.format(CSVFormat);
  file.close();
  
  // SOLUTION 
  file.open("mu1.csv");
  file << mu1.format(CSVFormat);
  file.close();

  file.open("mu2.csv");
  file << mu2.format(CSVFormat);
  file.close();

  file.open("mu3.csv");
  file << mu3.format(CSVFormat);
  file.close();

  file.open("mu4.csv");
  file << mu4.format(CSVFormat);
  file.close();

  file.open("mu5.csv");
  file << mu5.format(CSVFormat);
  file.close();
  
  // REFERENCE SOLUTION REDUCED
  file.open("mu_sub1.csv");
  file << mu_sub1.format(CSVFormat);
  file.close();

  file.open("mu_sub2.csv");
  file << mu_sub2.format(CSVFormat);
  file.close();

  file.open("mu_sub3.csv");
  file << mu_sub3.format(CSVFormat);
  file.close();

  file.open("mu_sub4.csv");
  file << mu_sub4.format(CSVFormat);
  file.close();

  file.open("mu_sub5.csv");
  file << mu_sub5.format(CSVFormat);
  file.close();

  return 0;
}
