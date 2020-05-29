/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 27.05.20
 *  @copyright ETH Zurich
 */

#include "norms.cc"
#include "strangsplitting.cc"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

#include <lf/io/io.h>

using namespace FisherKPP;

int main(int /*argc*/, char ** /*argv*/){
 
  Eigen::VectorXi numDofs(5); numDofs.setZero();
  Eigen::VectorXd meshsizes(5); meshsizes.setZero();
  Eigen::VectorXd eL2(5); eL2.setZero();
  Eigen::VectorXi m(5); m.setZero();
  Eigen::VectorXd tau(5); tau.setZero();

  Eigen::VectorXd sol1(39); sol1.setZero();
  Eigen::VectorXd sol2(125); sol2.setZero();
  Eigen::VectorXd sol3(444); sol3.setZero();
  Eigen::VectorXd sol4(1670); sol4.setZero();
  Eigen::VectorXd sol5(6474); sol5.setZero();

  Eigen::VectorXd cap1(39); cap1.setZero();
  cap1 = 0.8 * Eigen::VectorXd::Ones(39);
  
  Eigen::VectorXd cap2(125); cap2.setZero();
  cap2 = 0.8 * Eigen::VectorXd::Ones(125);

  Eigen::VectorXd cap3(444); cap3.setZero();
  cap3 = 0.8 * Eigen::VectorXd::Ones(444);

  Eigen::VectorXd cap4(1670); cap4.setZero();
  cap4 = 0.8 * Eigen::VectorXd::Ones(1670);

  Eigen::VectorXd cap5(6474); cap5.setZero();
  cap5 = 0.8 * Eigen::VectorXd::Ones(6474);

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
    std::cout << "N_dofs std" << N_dofs << std::endl;
	std::cout << "meshsize " << meshsizes << std::endl;

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
    double T = 1.; 			
    m(0) = 100;
    tau(0) = T/m(0);
	
	if(l > 1) {
	  tau(l-1) = (tau(0) / meshsizes(0)) * meshsizes(l-1);
	  m(l-1) = std::round(T/tau(l-1));
	}
    
	std::cout << "tau " << tau(l-1) << std::endl;
	std::cout << "m " << m(l-1) << std::endl;
    
	/* Now we may compute the solution */
    if(l == 1) {
	  StrangSplit StrangSplitter1(fe_space, T, m(0), lambda, c, h, L);
  	  sol1 = StrangSplitter1.Evolution(cap1, u0);
	} 
	else if(l == 2) {
  	  StrangSplit StrangSplitter2(fe_space, T, m(1), lambda, c, h, L);
  	  sol2 = StrangSplitter2.Evolution(cap2, u0);
	}
	else if(l == 3) {
	  StrangSplit StrangSplitter3(fe_space, T, m(2), lambda, c, h, L);
      sol3 = StrangSplitter3.Evolution(cap3, u0);
	}
	else if(l == 4) {
	  StrangSplit StrangSplitter4(fe_space, T, m(3), lambda, c, h, L);
  	  sol4 = StrangSplitter4.Evolution(cap4, u0);
	}
	else if(l == 5) {
  	  StrangSplit StrangSplitter5(fe_space, T, m(4), lambda, c, h, L);
      sol5 = StrangSplitter5.Evolution(cap5, u0);
	}

  }
  /* Compare the solution to the Reference solution on the finest mesh. */
  Eigen::VectorXd sol5_sub1(numDofs(0));
  Eigen::VectorXd sol5_sub2(numDofs(1));
  Eigen::VectorXd sol5_sub3(numDofs(2));
  Eigen::VectorXd sol5_sub4(numDofs(3));
  Eigen::VectorXd sol5_sub5(numDofs(4));

  sol5_sub1 = reduce(sol5, numDofs(0));
  sol5_sub2 = reduce(sol5, numDofs(1));
  sol5_sub3 = reduce(sol5, numDofs(2));
  sol5_sub4 = reduce(sol5, numDofs(3));
  sol5_sub5 = reduce(sol5, numDofs(4));

  eL2(0) = (sol1 - sol5_sub1).lpNorm<2>();
  eL2(1) = (sol2 - sol5_sub2).lpNorm<2>();
  eL2(2) = (sol3 - sol5_sub3).lpNorm<2>();
  eL2(3) = (sol4 - sol5_sub4).lpNorm<2>();
  eL2(4) = (sol5 - sol5_sub5).lpNorm<2>();

  std::cout << "errors computed" << std::endl;
 
  for(int l = 0; l < 5; l++) {

	std::cout << "numdofs for l = " << l << " : " << numDofs(l) << std::endl;
	std::cout << "meshsize for l = " << l << " : " << meshsizes(l) << std::endl;
    std::cout << "for l = " << l << "number of timesteps m : " << m(l) << std::endl;
	std::cout << "for l = " << l << "timestep size : " << tau(l) << std::endl;

    std::cout << "L2-error with respect to solution on finest mesh with smallest time step size for l = " << l << " : " << eL2(l) << std::endl;
  }

  // Define output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file;

  file.open("eL2.csv");
  file << eL2.format(CSVFormat);
  file.close();
  
  file.open("m.csv");
  file << m.format(CSVFormat);
  file.close();
  
  file.open("tau.csv");
  file << tau.format(CSVFormat);
  file.close();

  file.open("NumDofs.csv");
  file << numDofs.format(CSVFormat);
  file.close();

  file.open("meshsizes.csv");
  file << meshsizes.format(CSVFormat);
  file.close();
  
  /* SOLUTION */
  file.open("sol1.csv");
  file << sol1.format(CSVFormat);
  file.close();
  file.open("sol2.csv");
  file << sol2.format(CSVFormat);
  file.close();

  file.open("sol3.csv");
  file << sol3.format(CSVFormat);
  file.close();

  file.open("sol4.csv");
  file << sol4.format(CSVFormat);
  file.close();

  file.open("sol5.csv");
  file << sol5.format(CSVFormat);
  file.close();
  
  /* REFERENCE SOLUTION REDUCED */
  file.open("sol5_sub1.csv");
  file << sol5_sub1.format(CSVFormat);
  file.close();

  file.open("sol5_sub2.csv");
  file << sol5_sub2.format(CSVFormat);
  file.close();

  file.open("sol5_sub3.csv");
  file << sol5_sub3.format(CSVFormat);
  file.close();

  file.open("sol5_sub4.csv");
  file << sol5_sub4.format(CSVFormat);
  file.close();

  file.open("sol5_sub5.csv");
  file << sol5_sub5.format(CSVFormat);
  file.close();

  return 0;
}
