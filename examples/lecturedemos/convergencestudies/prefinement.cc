/**
 * @file prefinement.cc	
 * @brief Creates convergence plots for experiment 3.2.3.11
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include "felagrangeontria.h"
#include <iostream>




int main(int argc, char *argv[]) {
    for (unsigned p = 1 ; p < 5 ; ++p) {
	FeLagrangeONTria<double> refel(p);
	std::cout << refel.EvalReferenceShapeFunctions(refel.EvaluationNodes()) << std::endl;
    }
    return 0;
}
