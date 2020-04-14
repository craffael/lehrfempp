/**
 * @file prefinement.cc	
 * @brief Creates convergence plots for experiment 3.2.3.11
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include "felagrangeonquad.h"
#include "felagrangeontria.h"
#include "felagrangeonsegment.h"
#include <iostream>
#include <lf/uscalfe/uscalfe.h>




int main(int argc, char *argv[]) {
    const lf::uscalfe::FeLagrangeO1Quad<double> lf_o1;
    const lf::uscalfe::FeLagrangeO2Quad<double> lf_o2;
    const lf::uscalfe::FeLagrangeO3Quad<double> lf_o3;
    const FeLagrangeONQuad<double> own_o1(1);
    const FeLagrangeONQuad<double> own_o2(2);
    const FeLagrangeONQuad<double> own_o3(3);

    std::cout << "LF O1:" << std::endl;
    std::cout << lf_o1.EvalReferenceShapeFunctions(lf_o1.EvaluationNodes()) << std::endl;
    std::cout << "Own O1:" << std::endl;
    std::cout << own_o1.EvalReferenceShapeFunctions(own_o1.EvaluationNodes()) << std::endl << std::endl;
    std::cout << "LF O2:" << std::endl;
    std::cout << lf_o2.EvalReferenceShapeFunctions(lf_o2.EvaluationNodes()) << std::endl;
    std::cout << "Own O2:" << std::endl;
    std::cout << own_o2.EvalReferenceShapeFunctions(own_o2.EvaluationNodes()) << std::endl << std::endl;
    std::cout << "LF O3:" << std::endl;
    std::cout << lf_o3.EvalReferenceShapeFunctions(lf_o3.EvaluationNodes()) << std::endl;
    std::cout << "Own O3:" << std::endl;
    std::cout << own_o3.EvalReferenceShapeFunctions(own_o3.EvaluationNodes()) << std::endl << std::endl;

    return 0;
}
