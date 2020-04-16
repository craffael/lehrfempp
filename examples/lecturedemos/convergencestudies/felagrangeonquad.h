#ifndef EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONQUAD_H_
#define EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONQUAD_H_


#include <lf/uscalfe/uscalfe.h>
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <cmath>


template<typename SCALAR>
class FeLagrangeONQuad final : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
public:
    FeLagrangeONQuad(const FeLagrangeONQuad&) = default;
    FeLagrangeONQuad(FeLagrangeONQuad&&) = default;
    FeLagrangeONQuad& operator=(const FeLagrangeONQuad&) = default;
    FeLagrangeONQuad& operator=(FeLagrangeONQuad&&) = default;
    ~FeLagrangeONQuad() = default;

    FeLagrangeONQuad(unsigned degree) : degree_(degree), eval_nodes_(), ref_func_coeffs_() {
	eval_nodes_ = ComputeEvaluationNodes(degree);
	ref_func_coeffs_ = ComputePolyBasis(eval_nodes_).inverse().transpose();
    }

    [[nodiscard]] lf::base::RefEl RefEl() const override {
	return lf::base::RefEl::kQuad();
    }

    [[nodiscard]] unsigned Degree() const override {
	return degree_;
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
	return (degree_+1) * (degree_+1);
    }

    [[nodiscard]] lf::base::size_type NumRefShapeFunctions(lf::assemble::dim_t codim) const override {
	switch (codim) {
	    case 0:
		return (degree_-1) * (degree_-1);
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
	return ref_func_coeffs_ * poly_basis.transpose();
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
    Eigen::MatrixXd eval_nodes_;
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_func_coeffs_;

    Eigen::MatrixXd ComputeEvaluationNodes(unsigned p) const {
	Eigen::MatrixXd nodes(2, (p+1)*(p+1));
	// Add the evaluation nodes corresponding to the vertices
	nodes(0, 0) = 0;
	nodes(1, 0) = 0;
	nodes(0, 1) = 1;
	nodes(1, 1) = 0;
	nodes(0, 2) = 1;
	nodes(1, 2) = 1;
	nodes(0, 3) = 0;
	nodes(1, 3) = 1;
	// Add the evaluation nodes corresponding to the edges of the quad
	for (int i = 0 ; i < p-1 ; ++i) {
	    nodes(0, 4+i) = static_cast<double>(i+1) / p;
	    nodes(1, 4+i) = 0;
	}
	for (int i = 0 ; i < p-1 ; ++i) {
	    nodes(0, 3+p+i) = 1;
	    nodes(1, 3+p+i) = static_cast<double>(i+1) / p;
	}
	for (int i = 0 ; i < p-1 ; ++i) {
	    nodes(0, 2+p+p+i) = 1. - static_cast<double>(i+1) / p;
	    nodes(1, 2+p+p+i) = 1;
	}
	for (int i = 0 ; i < p-1 ; ++i) {
	    nodes(0, 1+p+p+p+i) = 0;
	    nodes(1, 1+p+p+p+i) = 1. - static_cast<double>(i+1) / p;
	}
	// Add the evaluation nodes corresponding to the interior of the quad
	for (int i = 0 ; i < p-1 ; ++i) {
	    for (int j = 0 ; j < p-1 ; ++j) {
		nodes(0, 4*p+(p-1)*i+j) = static_cast<double>(j+1) / p;
		nodes(1, 4*p+(p-1)*i+j) = static_cast<double>(i+1) / p;
	    }
	}
	return nodes;
    }

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasis(const Eigen::MatrixXd &refcoords) const {
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> poly(refcoords.cols(), (degree_+1)*(degree_+1));
	for (unsigned powx = 0 ; powx <= degree_ ; ++powx)  {
	    for (unsigned powy = 0 ; powy <= degree_ ; ++powy) {
		for (unsigned node_idx = 0 ; node_idx < refcoords.cols() ; ++node_idx) {
		    poly(node_idx, (degree_+1)*powx+powy) = std::pow(refcoords(0, node_idx), powx) * std::pow(refcoords(1, node_idx), powy);
		}
	    }
	}
	return poly;
    }

    std::pair<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>> ComputePolyBasisDerivative(const Eigen::MatrixXd &refcoords) const {
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dx(refcoords.cols(), (degree_+1)*(degree_+1));
	Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dy(refcoords.cols(), (degree_+1)*(degree_+1));
	for (unsigned powx = 0 ; powx <= degree_ ; ++powx) {
	    for (unsigned powy = 0 ; powy <= degree_ ; ++powy) {
		for (unsigned node_idx = 0 ; node_idx < refcoords.cols() ; ++node_idx) {
		    if (powx > 0) {
			coeffs_dx(node_idx, (degree_+1)*powx+powy) = powx * std::pow(refcoords(0, node_idx), powx-1) * std::pow(refcoords(1, node_idx), powy);
		    }
		    else {
			coeffs_dx(node_idx, (degree_+1)*powx+powy) = 0;
		    }
		    if (powy > 0) {
			coeffs_dy(node_idx, (degree_+1)*powx+powy) = powy * std::pow(refcoords(0, node_idx), powx) * std::pow(refcoords(1, node_idx), powy-1);
		    }
		    else {
			coeffs_dy(node_idx, (degree_+1)*powx+powy) = 0;
		    }
		}
	    }
	}
	return {coeffs_dx, coeffs_dy};
    }
};


#endif // EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FELAGRANGEONQUAD_H_
