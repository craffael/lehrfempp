/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 01.04.20
 *  @copyright ETH Zurich
 */

#include <memory>
#include <utility>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace FisherKPP {

class StrangSplit {

  public:
    /* Disabled constructors */
	StrangSplit() = delete;
    StrangSplit(const StrangSplit &) = delete;
    StrangSplit(StrangSplit &&) = delete;
    StrangSplit &operator=(const StrangSplit &) = delete;
    StrangSplit &operator=(const StrangSplit &&) = delete;
	/* Main constructor */
	template <typename DIFF_COEFF, typename NONLOC_BC>
	explicit StrangSplit(const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space,
						 double T, unsigned int m, double lambda, DIFF_COEFF &&c, NONLOC_BC &&h, Eigen::MatrixXd L);
	/* Destructor */
	virtual ~StrangSplit() = default;

	/* Member Functions */
	Eigen::VectorXd diffusionEvolutionOperator(double tau, const Eigen::VectorXd &mu);
	
	//template <typename CAR_CAP>
	Eigen::VectorXd Evolution(Eigen::VectorXd &cap, const Eigen::VectorXd &mu);

  private:
    /* Finite Element Space */
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
	/* Final Time */
	double T_;
	/* Number of Timesteps */
    unsigned int m_;
	/* Growth Factor */
    double lambda_;
	/* coefficient for SDIRK-2 Butcher Tableau */
    double kappa_;
	/* Galerkin matrix corresponding to the negative Laplacian with Robin Boundary Conditions */
    Eigen::SparseMatrix<double> A_;
    /* Galerkin matrix for the Mass Matrix */
    Eigen::SparseMatrix<double> M_;
    /* Precompute LU decomposition needed for time stepping */
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

};

} /* namespace FisherKPP */
