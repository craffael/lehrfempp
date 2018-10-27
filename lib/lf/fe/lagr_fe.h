#ifndef LF_FE
#define LF_FE
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structures representing simple finite elements
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "dofhandler.h"

namespace lf::assemble {

  template<typename SCALAR,base::RefElType REF_EL, int Degree>
  class 2DLagrangeFiniteElement {
  public:
    constexpr size_type NumRefShapeFunctions() const { return (p+1)*(p+2)/2; }
    
    std::array<Eigen::Matrix<SCALAR,Eigen::Dynamic,1>,NumRefShapeFunctions()>
    EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords);

    std::array<Eigen::Matrix<SCALAR,2,Eigen::Dynamic>,NumRefShapeFunctions()>
    GradientsReferenceShapeFunctions(const Eigen::MatrixXd &refcoords);

    constexpr UniformFEDofHandler::dof_map_t RefShapeFunctionPattern() const {
      switch(p) {
      case 1: { 
	return (UniformFEDofHandler::dof_map_t{{lf::base::RefEl::kPoint(),1}});
	break;
      }
      case 2: {
	return (UniformFEDofHandler::dof_map_t{{lf::base::RefEl::kPoint(),1},{lf::base::RefEl::kSegment(),1});
	break;
      }
      

    
  }
}

#endif
