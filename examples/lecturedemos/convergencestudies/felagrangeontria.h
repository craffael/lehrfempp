#ifndef EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONTRIA_H_
#define EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONTRIA_H_



#include <lf/uscalfe/uscalfe.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <cmath>


template<typename SCALAR>
class FeLagrangeONTria final : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
public:
    FeLagrangeONTria(const FeLagrangeONTria&) = default;
    FeLagrangeONTria(FeLagrangeONTria&&) = default;
    FeLagrangeONTria& operator=(const FeLagrangeONTria&) = default;
    FeLagrangeONTria& operator=(FeLagrangeONTria&&) = default;
    ~FeLagrangeONTria() = default;

    FeLagrangeONTria(unsigned degree) : degree_(degree), eval_nodes_(), ref_func_coeffs_() {
	eval_nodes_ = ComputeEvaluationNodes(degree);
	ref_func_coeffs_ = ComputePolyBasis(eval_nodes_).inverse().transpose();
    }

    [[nodiscard]] lf::base::RefEl RefEl() const override {
	return lf::base::RefEl::kTria();
    }

    [[nodiscard]] unsigned Degree() const override {
	return degree_;
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
	return (degree_+1) * (degree_+2) / 2;
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions(lf::assemble::dim_t codim) const override {
	switch (codim) {
	    case 0:
		if (degree_ <= 2) {
		    return 0;
		}
		else {
		    return (degree_-2) * (degree_-1) / 2;
		}
	    case 1:
		return degree_ - 1;
	    case 2:
		return 1;
	    default:
		LF_ASSERT_MSG(false, "Illegal codim " << codim);
	}
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions(lf::assemble::dim_t codim, lf::base::sub_idx_t /*subidx*/) const override {
	return NumRefShapeFunctions(codim);
    }

    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
	const auto poly_basis = ComputePolyBasis(refcoords);
	return ref_func_coeffs_ * poly_basis;
    }
    
    [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> GradientsReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> grads(NumRefShapeFunctions(), 2*refcoords.cols());
	const auto [basis_dx, basis_dy] = ComputePolyBasisDerivative(refcoords);
	const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dx = ref_func_coeffs_ * basis_dx.transpose();
	const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dy = ref_func_coeffs_ * basis_dy.transpose();
	for (int refcoord_idx = 0 ; refcoord_idx < refcoords.cols() ; ++refcoord_idx) {
	    grads.col(2*refcoord_idx+0) = dx.col(refcoord_idx);
	    grads.col(2*refcoord_idx+1) = dy.col(refcoord_idx);
	}
	return grads;
    }

    [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
	return eval_nodes_;
    }

    [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
	return NumRefShapeFunctions();
    }


private:
    unsigned degree_;
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> eval_nodes_;
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_func_coeffs_;

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputeEvaluationNodes(unsigned p) const {
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> eval_nodes(2, (p+1)*(p+2)/2);
	// Add the evaluation nodes corresponding to the vertices of the triangle
	eval_nodes(0, 0) = 0;
	eval_nodes(1, 0) = 0;
	eval_nodes(0, 1) = 1;
	eval_nodes(1, 1) = 0;
	eval_nodes(0, 2) = 0;
	eval_nodes(1, 2) = 1;
	// Add the evaluation nodes corresponding to the edges of the triangle
	for (int i = 0 ; i < p-1 ; ++i) {
	    eval_nodes(0, 3+i) = static_cast<SCALAR>(i+1) / p;
	    eval_nodes(1, 3+i) = 0;
	}
	for (int i = 0 ; i < p-1 ; ++i) {
	    eval_nodes(0, 2+p+i) = 1. - (static_cast<SCALAR>(i+1) / p);
	    eval_nodes(1, 2+p+i) = static_cast<SCALAR>(i+1) / p;
	}
	for (int i = 0 ; i < p-1 ; ++i) {
	    eval_nodes(0, 1+p+p+i) = 0;
	    eval_nodes(1, 1+p+p+i) = 1. - (static_cast<SCALAR>(i+1) / p);
	}
	// Add the evaluation nodes corresponding to the interior of the triangle
	if (p > 2) {
	    for (int i = 0 ; i < p-2 ; ++i) {
		for (int j = 0 ; j <= i ; ++j) {
		    eval_nodes(0, (3*p)+(i*(i+1)/2)+j) = static_cast<SCALAR>(i-j+1) / p;
		    eval_nodes(1, (3*p)+(i*(i+1)/2)+j) = static_cast<SCALAR>(j+1) / p;
		}
	    }
	}
	return eval_nodes;
    }

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasis(const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> &refcoords) const {
	// The coefficients are ordered x^0y^0, x^0y^1, ..., x^0y^p, x^1y^0, x^1y^1, ..., x^1, y^(p-1), ... x^(p-1)y^0, x^(p-1)y^1, x^py^0
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs(refcoords.cols(), (degree_+1)*(degree_+2)/2);
	unsigned coeff_idx = 0;
	for (unsigned powx = 0 ; powx <= degree_ ; ++powx) {
	    for (unsigned powy = 0 ; powy <= degree_-powx ; ++powy) {
		for (unsigned node_idx = 0 ; node_idx < refcoords.cols() ; ++node_idx) {
		    coeffs(node_idx, coeff_idx) = std::pow(refcoords(0, node_idx), powx) * std::pow(refcoords(1, node_idx), powy);
		}
		++coeff_idx;
	    }
	}
	return coeffs;
    }

    std::pair<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>> ComputePolyBasisDerivative(const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> &refcoords) const {
	// The coefficients are ordered x^0y^0, x^0y^1, ..., x^0y^p, x^1y^0, x^1y^1, ..., x^1, y^(p-1), ... x^(p-1)y^0, x^(p-1)y^1, x^py^0
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dx(refcoords.cols(), (degree_+1)*(degree_+2)/2);
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dy(refcoords.cols(), (degree_+1)*(degree_+2)/2);
	unsigned coeff_idx = 0;
	for (unsigned powx = 0 ; powx <= degree_ ; ++powx) {
	    for (unsigned powy = 0 ; powy <= degree_-powx ; ++powy) {
		for (unsigned node_idx = 0 ; node_idx < refcoords.cols() ; ++node_idx) {
		    if (powx > 0) {
			coeffs_dx(node_idx, coeff_idx) = powx * std::pow(refcoords(0, node_idx), powx-1) * std::pow(refcoords(1, node_idx), powy);
		    }
		    else {
			coeffs_dx(node_idx, coeff_idx) = 0;
		    }
		    if (powy > 0) {
			coeffs_dy(node_idx, coeff_idx) = powy * std::pow(refcoords(0, node_idx), powx) * std::pow(refcoords(1, node_idx), powy-1);
		    }
		    else {
			coeffs_dy(node_idx, coeff_idx) = 0;
		    }
		}
		++coeff_idx;
	    }
	}
	return {coeffs_dx, coeffs_dy};
    }
};



#endif // EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONTRIA_H_
